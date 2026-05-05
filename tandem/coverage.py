"""Coverage analysis for Tandem module 2.

Coverage detection: simple threshold-based with two-pass merging.
Coverage plot: interactive HTML/JS with zoom, pan, and coordinate readout.
Users inspect the plot visually and provide coordinates to the junction
discovery pipeline.
"""

import logging
import os

import numpy as np

from . import utils

logger = logging.getLogger("tandem")


# =============================================================================
# Coverage calculation
# =============================================================================

def calculate_coverage_from_bam(bam_path, seq_id=None):
    """Calculate per-position coverage from a BAM file using samtools depth."""
    cmd = ["samtools", "depth", "-a", "-J"]
    if seq_id:
        cmd += ["-r", seq_id]
    cmd.append(str(bam_path))

    result = utils.run_command(cmd, description="Calculating coverage")

    coverage = {}
    for line in result.stdout.split("\n"):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        sid = parts[0]
        pos = int(parts[1]) - 1
        depth = int(parts[2])

        if sid not in coverage:
            coverage[sid] = []
        while len(coverage[sid]) <= pos:
            coverage[sid].append(0)
        coverage[sid][pos] = depth

    for sid in coverage:
        coverage[sid] = np.array(coverage[sid], dtype=np.int32)

    return coverage


def sliding_window_coverage(coverage_array, window_size=200):
    """Calculate mean coverage in sliding windows."""
    n = len(coverage_array)
    if n < window_size:
        return np.array([np.mean(coverage_array)]), np.array([0])

    n_windows = n - window_size + 1
    cumsum = np.cumsum(coverage_array, dtype=np.float64)
    cumsum = np.insert(cumsum, 0, 0)
    windowed = (cumsum[window_size:] - cumsum[:-window_size]) / window_size
    starts = np.arange(n_windows)

    return windowed, starts


# =============================================================================
# Threshold-based detection with two-pass merging
# =============================================================================

def detect_elevated_regions(coverage_array, window_size=200,
                            min_fold_change=1.7, min_region_size=500,
                            merge_fold_threshold=1.3,
                            merge_distance=3000):
    """Detect regions with elevated coverage using threshold approach.

    Returns:
        (list of region dicts, genome_median)
    """
    windowed, starts = sliding_window_coverage(coverage_array, window_size)

    genome_median = float(np.median(windowed))
    genome_mad = float(np.median(np.abs(windowed - genome_median)))

    if genome_median == 0:
        logger.warning("Genome median coverage is 0. No reads mapped?")
        return [], 0.0

    threshold = genome_median * min_fold_change
    logger.info(
        f"  Coverage: median={genome_median:.1f}x, "
        f"MAD={genome_mad:.1f}, threshold={threshold:.1f}x "
        f"({min_fold_change:.1f}-fold)"
    )

    elevated_mask = windowed >= threshold
    raw_regions = _mask_to_regions(elevated_mask, starts, window_size)
    raw_regions = [r for r in raw_regions if r[1] - r[0] >= min_region_size]

    logger.info(f"  Found {len(raw_regions)} raw elevated regions before merging")

    merge_threshold = genome_median * merge_fold_threshold
    merged = _coverage_aware_merge(raw_regions, coverage_array,
                                    merge_threshold, window_size)

    logger.info(
        f"  After coverage-aware merging (gap threshold {merge_fold_threshold:.1f}x = "
        f"{merge_threshold:.1f}x): {len(merged)} regions"
    )

    merged = _merge_regions(merged, merge_distance)

    logger.info(
        f"  After distance merging ({merge_distance:,} bp): {len(merged)} regions"
    )

    results = []
    for start, end in merged:
        start = int(start)
        end = int(min(end, len(coverage_array)))
        region_cov = coverage_array[start:end]
        mean_cov = float(np.mean(region_cov))
        fc = mean_cov / genome_median if genome_median > 0 else 0

        results.append({
            "start": start,
            "end": end,
            "mean_coverage": round(mean_cov, 1),
            "fold_change": round(fc, 2),
            "size": end - start,
            "genome_median_coverage": round(genome_median, 1),
        })

    logger.info(f"  Detected {len(results)} elevated regions")
    return results, genome_median


# =============================================================================
# Helper functions
# =============================================================================

def _mask_to_regions(mask, starts, window_size):
    """Convert boolean mask to list of (start, end) regions."""
    regions = []
    in_region = False
    region_start = 0

    for i, is_elevated in enumerate(mask):
        if is_elevated and not in_region:
            region_start = starts[i]
            in_region = True
        elif not is_elevated and in_region:
            region_end = starts[i] + window_size
            regions.append((region_start, region_end))
            in_region = False

    if in_region:
        region_end = starts[len(mask) - 1] + window_size
        regions.append((region_start, region_end))

    return regions


def _coverage_aware_merge(regions, coverage_array, merge_threshold,
                           window_size=200):
    """Merge adjacent regions if coverage in the gap stays elevated."""
    if len(regions) <= 1:
        return regions

    sorted_regions = sorted(regions, key=lambda r: r[0])

    changed = True
    while changed:
        changed = False
        new_regions = [sorted_regions[0]]

        for start, end in sorted_regions[1:]:
            prev_start, prev_end = new_regions[-1]
            gap_start = int(prev_end)
            gap_end = int(start)

            if gap_end <= gap_start:
                new_regions[-1] = (prev_start, max(prev_end, end))
                changed = True
                continue

            gap_cov = coverage_array[gap_start:gap_end]
            gap_mean = float(np.mean(gap_cov)) if len(gap_cov) > 0 else 0.0

            if gap_mean >= merge_threshold:
                new_regions[-1] = (prev_start, max(prev_end, end))
                changed = True
            else:
                new_regions.append((start, end))

        sorted_regions = new_regions

    return sorted_regions


def _merge_regions(regions, max_gap):
    """Merge regions that are within max_gap of each other."""
    if not regions:
        return []

    sorted_regions = sorted(regions, key=lambda r: r[0])
    merged = [sorted_regions[0]]

    for start, end in sorted_regions[1:]:
        prev_start, prev_end = merged[-1]
        if start - prev_end <= max_gap:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))

    return merged


# =============================================================================
# Interactive coverage plot (self-contained HTML/JS, no CDN)
# =============================================================================

def generate_coverage_plot(coverage_array, seq_id, output_dir,
                           window_size=200, genome_median=None,
                           elevated_regions=None, min_fold_change=1.7):
    """Generate an interactive HTML coverage plot.

    Features:
    - Mouse wheel zoom (centered on cursor)
    - Click-drag to pan
    - Crosshair with exact position and coverage readout
    - Threshold and median lines
    - Highlighted elevated regions
    - Fully self-contained (no CDN, works offline on HPC)

    Uses 200bp windows for high-resolution coordinate readout when
    zoomed in. Users can identify boundaries to ±200bp accuracy,
    matching the default -flank 200 for junction discovery.
    """
    windowed, starts = sliding_window_coverage(coverage_array, window_size)

    if genome_median is None:
        genome_median = float(np.median(windowed))

    threshold = genome_median * min_fold_change
    genome_len = len(coverage_array)

    # Keep all data points — 200bp windows on 6.7Mb = ~33k points
    # Canvas handles this fine; enables precise coordinate readout when zoomed
    max_points = 50000
    if len(windowed) > max_points:
        step = len(windowed) // max_points
        plot_pos = starts[::step].tolist()
        plot_cov = [round(float(c), 1) for c in windowed[::step]]
    else:
        plot_pos = starts.tolist()
        plot_cov = [round(float(c), 1) for c in windowed]

    # Regions as JS array
    regions_js = "[]"
    if elevated_regions:
        items = []
        for r in elevated_regions:
            items.append(
                f'{{s:{r["start"]},e:{r["end"]},fc:{r["fold_change"]:.2f},'
                f'sz:{r["size"]}}}'
            )
        regions_js = "[" + ",".join(items) + "]"

    # Regions table HTML
    rtable = ""
    if elevated_regions:
        rows = "".join(
            f'<tr onclick="zoomTo({r["start"]},{r["end"]})" style="cursor:pointer">'
            f'<td>{i+1}</td><td>{r["start"]:,}</td><td>{r["end"]:,}</td>'
            f'<td>{r["size"]:,}</td><td>{r["fold_change"]:.2f}x</td></tr>\n'
            for i, r in enumerate(elevated_regions)
        )
        rtable = f"""
    <h3>Detected Regions (click to zoom)</h3>
    <table class="rt"><thead><tr><th>#</th><th>Start</th><th>End</th>
    <th>Size (bp)</th><th>Fold Change</th></tr></thead>
    <tbody>{rows}</tbody></table>"""

    html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8">
<title>Tandem Coverage - {seq_id}</title>
<style>
*{{margin:0;padding:0;box-sizing:border-box}}
body{{font-family:-apple-system,BlinkMacSystemFont,sans-serif;background:#f5f5f5;user-select:none}}
.hd{{background:#2c3e50;color:#fff;padding:14px 24px}}
.hd h1{{font-size:20px;display:inline}}.hd .s{{font-size:13px;opacity:.8;margin-left:12px}}
.c{{max-width:1600px;margin:12px auto;padding:0 16px}}
.info{{background:#fff;padding:10px 16px;border-radius:6px;margin-bottom:10px;font-size:13px;
box-shadow:0 1px 3px rgba(0,0,0,.1);display:flex;gap:20px;flex-wrap:wrap;align-items:center}}
.info label{{font-weight:600;color:#666}}
.cb{{background:#fff;border-radius:6px;padding:8px;box-shadow:0 1px 3px rgba(0,0,0,.1);position:relative}}
canvas{{display:block;cursor:crosshair}}
.coords{{position:absolute;top:8px;right:12px;background:rgba(255,255,255,0.92);
padding:4px 10px;border-radius:4px;font-size:12px;font-family:monospace;
box-shadow:0 1px 3px rgba(0,0,0,.15);pointer-events:none}}
.controls{{display:flex;gap:8px;margin:8px 0;font-size:13px;align-items:center}}
.controls button{{padding:4px 12px;border:1px solid #ccc;border-radius:4px;background:#fff;cursor:pointer}}
.controls button:hover{{background:#ecf0f1}}
.tip{{background:#eef6ff;border-left:4px solid #4A90D9;padding:10px 16px;margin-top:10px;
font-size:13px;border-radius:0 4px 4px 0}}
.tip code{{background:#ddd;padding:1px 4px;border-radius:3px;font-size:12px}}
h3{{margin-top:14px;font-size:14px;margin-bottom:6px}}
.rt{{width:100%;border-collapse:collapse;font-size:12px;background:#fff;border-radius:6px;
box-shadow:0 1px 3px rgba(0,0,0,.1)}}
.rt th{{background:#ecf0f1;padding:6px 10px;text-align:left;font-weight:600;color:#555}}
.rt td{{padding:5px 10px;border-top:1px solid #eee}}
.rt tr:hover{{background:#f0f7ff}}
</style></head><body>
<div class="hd"><h1>Tandem</h1><span class="s">{seq_id} &middot; {genome_len:,} bp</span></div>
<div class="c">
<div class="info">
<span><label>Median:</label> {genome_median:.1f}x</span>
<span><label>Threshold:</label> {threshold:.1f}x ({min_fold_change:.1f}&times;)</span>
<span><label>Window:</label> {window_size:,} bp</span>
<span><label>Regions:</label> {len(elevated_regions) if elevated_regions else 0}</span>
</div>
<div class="controls">
<button onclick="resetZoom()">Reset Zoom</button>
<button onclick="zoomIn()">Zoom In</button>
<button onclick="zoomOut()">Zoom Out</button>
<span id="rangeLabel" style="color:#888;margin-left:8px"></span>
</div>
<div class="cb">
<canvas id="cv" width="1500" height="450"></canvas>
<div class="coords" id="coords">Position: — &nbsp; Coverage: —</div>
</div>
<div class="tip">
<strong>Scroll</strong> to zoom, <strong>drag</strong> to pan, <strong>click region table</strong> to zoom to region.<br>
Then run: <code>tandem -r ref.fna -i R1.fq -I R2.fq -iso -s START -e END -flank 200 -t THREADS -o output</code>
</div>
{rtable}
</div>
<script>
const P={plot_pos};
const C={plot_cov};
const MED={genome_median:.1f};
const THR={threshold:.1f};
const REGIONS={regions_js};
const GLEN={genome_len};

const cv=document.getElementById('cv');
const ctx=cv.getContext('2d');
const W=cv.width, H=cv.height;
const ML=65,MR=25,MT=25,MB=45;
const PW=W-ML-MR, PH=H-MT-MB;

let viewStart=P[0], viewEnd=P[P.length-1];
let maxCov=Math.max(...C)*1.15;

function sx(p){{ return ML+(p-viewStart)/(viewEnd-viewStart)*PW; }}
function sy(c){{ return MT+PH-c/maxCov*PH; }}
function invX(px){{ return viewStart+(px-ML)/PW*(viewEnd-viewStart); }}
function invY(py){{ return (MT+PH-py)/PH*maxCov; }}

function draw(){{
  ctx.clearRect(0,0,W,H);

  // Regions
  ctx.fillStyle='rgba(74,144,217,0.12)';
  ctx.strokeStyle='rgba(74,144,217,0.35)';
  ctx.lineWidth=0.5;
  for(const r of REGIONS){{
    const x1=sx(r.s),x2=sx(r.e);
    if(x2<ML||x1>ML+PW) continue;
    const rx=Math.max(ML,x1), rw=Math.min(ML+PW,x2)-rx;
    ctx.fillRect(rx,MT,rw,PH);
    ctx.strokeRect(rx,MT,rw,PH);
  }}

  // Grid
  ctx.strokeStyle='#eee';ctx.lineWidth=0.5;
  for(let i=0;i<=5;i++){{
    const y=sy(maxCov*i/5);
    ctx.beginPath();ctx.moveTo(ML,y);ctx.lineTo(ML+PW,y);ctx.stroke();
  }}

  // Coverage line
  ctx.beginPath();
  ctx.strokeStyle='#2c3e50';ctx.lineWidth=1;
  let started=false;
  for(let i=0;i<P.length;i++){{
    const x=sx(P[i]),y=sy(C[i]);
    if(x<ML-2||x>ML+PW+2) continue;
    if(!started){{ctx.moveTo(x,y);started=true;}} else ctx.lineTo(x,y);
  }}
  ctx.stroke();

  // Fill
  ctx.beginPath();
  started=false;
  let firstX=0,lastX=0;
  for(let i=0;i<P.length;i++){{
    const x=sx(P[i]),y=sy(C[i]);
    if(x<ML-2||x>ML+PW+2) continue;
    if(!started){{ctx.moveTo(x,y);firstX=x;started=true;}} else ctx.lineTo(x,y);
    lastX=x;
  }}
  ctx.lineTo(lastX,sy(0));ctx.lineTo(firstX,sy(0));ctx.closePath();
  ctx.fillStyle='rgba(44,62,80,0.06)';ctx.fill();

  // Median line
  const my=sy(MED);
  ctx.strokeStyle='#27ae60';ctx.lineWidth=1;ctx.setLineDash([5,4]);
  ctx.beginPath();ctx.moveTo(ML,my);ctx.lineTo(ML+PW,my);ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle='#27ae60';ctx.font='10px sans-serif';
  ctx.fillText('Median '+MED.toFixed(1)+'x',ML+3,my-3);

  // Threshold line
  const ty=sy(THR);
  ctx.strokeStyle='#e74c3c';ctx.lineWidth=1;ctx.setLineDash([5,4]);
  ctx.beginPath();ctx.moveTo(ML,ty);ctx.lineTo(ML+PW,ty);ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle='#e74c3c';ctx.font='10px sans-serif';
  ctx.fillText('Threshold '+THR.toFixed(1)+'x',ML+3,ty-3);

  // Axes
  ctx.strokeStyle='#333';ctx.lineWidth=1;
  ctx.beginPath();ctx.moveTo(ML,MT);ctx.lineTo(ML,MT+PH);ctx.lineTo(ML+PW,MT+PH);ctx.stroke();

  // X ticks
  ctx.fillStyle='#666';ctx.font='11px sans-serif';ctx.textAlign='center';
  const range=viewEnd-viewStart;
  const nxt=Math.min(15, Math.max(5, Math.floor(PW/80)));
  const step=range/nxt;
  for(let i=0;i<=nxt;i++){{
    const v=viewStart+i*step;
    const x=sx(v);
    ctx.beginPath();ctx.moveTo(x,MT+PH);ctx.lineTo(x,MT+PH+4);ctx.stroke();
    let lb;
    if(range>1e6) lb=(v/1e6).toFixed(2)+'M';
    else if(range>1e4) lb=(v/1e3).toFixed(0)+'k';
    else lb=v.toFixed(0);
    ctx.fillText(lb,x,MT+PH+16);
  }}

  // Y ticks
  ctx.textAlign='right';
  for(let i=0;i<=5;i++){{
    const cv=maxCov*i/5;
    const y=sy(cv);
    ctx.beginPath();ctx.moveTo(ML-4,y);ctx.lineTo(ML,y);ctx.stroke();
    ctx.fillText(cv.toFixed(0)+'x',ML-6,y+4);
  }}

  // Axis labels
  ctx.fillStyle='#666';ctx.font='12px sans-serif';ctx.textAlign='center';
  ctx.fillText('Genome position (bp)',ML+PW/2,H-5);
  ctx.save();ctx.translate(12,MT+PH/2);ctx.rotate(-Math.PI/2);
  ctx.fillText('Coverage (x)',0,0);ctx.restore();

  // Range label
  document.getElementById('rangeLabel').textContent=
    'Viewing: '+(viewStart/1e6).toFixed(3)+'M - '+(viewEnd/1e6).toFixed(3)+'M ('+
    ((viewEnd-viewStart)/1e3).toFixed(1)+'kb)';
}}

// Mouse interaction
let isDragging=false, dragStartX=0, dragViewStart=0, dragViewEnd=0;

cv.addEventListener('wheel',function(e){{
  e.preventDefault();
  const rect=cv.getBoundingClientRect();
  const mx=e.clientX-rect.left;
  const pos=invX(mx);
  const factor=e.deltaY>0?1.3:1/1.3;
  const newRange=(viewEnd-viewStart)*factor;
  const frac=(mx-ML)/PW;
  viewStart=pos-frac*newRange;
  viewEnd=pos+(1-frac)*newRange;
  viewStart=Math.max(0,viewStart);
  viewEnd=Math.min(GLEN,viewEnd);
  if(viewEnd-viewStart<500){{viewEnd=viewStart+500;}}
  draw();
}});

cv.addEventListener('mousedown',function(e){{
  isDragging=true;
  dragStartX=e.clientX;
  dragViewStart=viewStart;
  dragViewEnd=viewEnd;
  cv.style.cursor='grabbing';
}});

window.addEventListener('mousemove',function(e){{
  if(isDragging){{
    const dx=e.clientX-dragStartX;
    const range=dragViewEnd-dragViewStart;
    const shift=-dx/PW*range;
    viewStart=Math.max(0,dragViewStart+shift);
    viewEnd=Math.min(GLEN,dragViewEnd+shift);
    draw();
  }}
  // Crosshair coordinates
  const rect=cv.getBoundingClientRect();
  const mx=e.clientX-rect.left;
  const my=e.clientY-rect.top;
  if(mx>=ML&&mx<=ML+PW&&my>=MT&&my<=MT+PH){{
    const pos=invX(mx);
    const cov=invY(my);
    // Find nearest data point for actual coverage
    let nearest=0,bestDist=Infinity;
    for(let i=0;i<P.length;i++){{
      const d=Math.abs(P[i]-pos);
      if(d<bestDist){{bestDist=d;nearest=i;}}
    }}
    const actualCov=C[nearest];
    document.getElementById('coords').innerHTML=
      'Position: <b>'+Math.round(pos).toLocaleString()+'</b> bp &nbsp; '+
      'Coverage: <b>'+actualCov.toFixed(1)+'x</b>';
  }}
}});

window.addEventListener('mouseup',function(){{
  isDragging=false;
  cv.style.cursor='crosshair';
}});

function resetZoom(){{viewStart=P[0];viewEnd=P[P.length-1];draw();}}
function zoomIn(){{
  const mid=(viewStart+viewEnd)/2;
  const range=(viewEnd-viewStart)/3;
  viewStart=Math.max(0,mid-range);viewEnd=Math.min(GLEN,mid+range);draw();
}}
function zoomOut(){{
  const mid=(viewStart+viewEnd)/2;
  const range=(viewEnd-viewStart)*1.5;
  viewStart=Math.max(0,mid-range);viewEnd=Math.min(GLEN,mid+range);draw();
}}
function zoomTo(s,e){{
  const pad=(e-s)*0.2;
  viewStart=Math.max(0,s-pad);viewEnd=Math.min(GLEN,e+pad);draw();
}}

draw();
</script>
</body></html>"""

    output_path = os.path.join(output_dir, f"coverage_plot_{seq_id}.html")
    with open(output_path, "w") as f:
        f.write(html)

    logger.info(f"  Coverage plot saved to {output_path}")
    return output_path
