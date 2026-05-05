"""HTML report generator for Tandem.

Generates an interactive HTML report inspired by breseq output format:
- Summary table of all detected duplications/junctions
- Click-through to detailed junction alignment views
- Mechanism classification coloring
"""

import json
import logging
import os
from pathlib import Path

logger = logging.getLogger("tandem")


# =============================================================================
# Color scheme for mechanisms
# =============================================================================

MECHANISM_COLORS = {
    "NHEJ": "#4A90D9",
    "MMEJ": "#F5A623",
    "SSA": "#D0021B",
    "HR": "#7ED321",
    "unknown": "#999999",
}

CONFIDENCE_BADGES = {
    "high": '<span class="badge badge-high">high</span>',
    "medium": '<span class="badge badge-medium">medium</span>',
    "low": '<span class="badge badge-low">low</span>',
}


# =============================================================================
# Report generation
# =============================================================================

def generate_report(results, module_name, output_dir, ref_name="",
                    extra_info=None):
    """Generate an HTML report for Tandem results.

    Args:
        results: list of result dicts (from any module)
        module_name: 'module1', 'module2', or 'module3'
        output_dir: output directory
        ref_name: reference genome name for display
        extra_info: dict of additional metadata to display
    """
    output_dir = Path(output_dir)
    report_path = output_dir / "tandem_report.html"

    if module_name == "module1":
        title = "Tandem - Reference Genome Analysis"
        table_html = _build_module1_table(results)
        detail_html = _build_module1_details(results)
    elif module_name == "module2":
        title = "Tandem - Isolate Junction Discovery"
        table_html = _build_module2_table(results)
        detail_html = _build_module2_details(results)
    elif module_name == "module3":
        title = "Tandem - Population Junction Quantification"
        table_html = _build_module3_table(results)
        detail_html = _build_module3_details(results)
    else:
        title = "Tandem - Results"
        table_html = "<p>No results to display.</p>"
        detail_html = ""

    # Build full HTML
    html = _build_html_page(title, ref_name, table_html, detail_html,
                            module_name, extra_info)

    with open(report_path, "w") as f:
        f.write(html)

    logger.info(f"  HTML report saved to {report_path}")
    return str(report_path)


# =============================================================================
# Module-specific tables
# =============================================================================

def _build_module1_table(results):
    """Build summary table for module 1 results."""
    if not results:
        return "<p>No tandem duplications detected.</p>"

    rows = []
    for i, r in enumerate(results):
        is_hr = r.get("is_hr_signature", False)
        mh_bp = r.get("microhomology_bp", 0)
        hr_label = "HR" if is_hr else "non-HR"
        color = MECHANISM_COLORS.get("HR" if is_hr else "NHEJ", "#999")

        row = f"""
        <tr class="clickable-row" data-target="detail-{i}">
            <td>{i + 1}</td>
            <td>{r.get('seq_id', '')}</td>
            <td>{r.get('copy1_start', ''):,}</td>
            <td>{r.get('copy1_end', ''):,}</td>
            <td>{r.get('copy2_start', ''):,}</td>
            <td>{r.get('copy2_end', ''):,}</td>
            <td>{r.get('size', 0):,} bp</td>
            <td>{r.get('distance', 0):,} bp</td>
            <td>{r.get('identity', 0):.1f}%</td>
            <td><span class="mechanism-tag" style="background:{color}">{hr_label}</span></td>
            <td>{mh_bp} bp</td>
        </tr>"""
        rows.append(row)

    return f"""
    <table class="results-table" id="summary-table">
        <thead>
            <tr>
                <th>#</th>
                <th>Sequence</th>
                <th>Copy1 Start</th>
                <th>Copy1 End</th>
                <th>Copy2 Start</th>
                <th>Copy2 End</th>
                <th>Size</th>
                <th>Distance</th>
                <th>Identity</th>
                <th>HR Status</th>
                <th>Microhomology</th>
            </tr>
        </thead>
        <tbody>
            {''.join(rows)}
        </tbody>
    </table>
    <p class="table-note">Click any row for detailed junction view.</p>
    """


def _build_module2_table(results):
    """Build summary table for module 2 results."""
    if not results:
        return "<p>No junctions confirmed.</p>"

    rows = []
    for i, r in enumerate(results):
        is_hr = r.get("is_hr_signature", False)
        mh_bp = r.get("microhomology_bp", 0)
        hr_label = "HR" if is_hr else "non-HR"
        color = MECHANISM_COLORS.get("HR" if is_hr else "NHEJ", "#999")

        row = f"""
        <tr class="clickable-row" data-target="detail-{i}">
            <td>{i + 1}</td>
            <td>{r.get('seq_id', '')}</td>
            <td>{r.get('dup_start', 0):,}</td>
            <td>{r.get('dup_end', 0):,}</td>
            <td>{r.get('dup_size', 0):,} bp</td>
            <td><span class="mechanism-tag" style="background:{color}">{hr_label}</span></td>
            <td>{mh_bp} bp</td>
            <td>{r.get('spanning_reads', 0)}</td>
            <td>{r.get('hq_reads', 0)}</td>
        </tr>"""
        rows.append(row)

    return f"""
    <table class="results-table" id="summary-table">
        <thead>
            <tr>
                <th>#</th>
                <th>Sequence</th>
                <th>Dup Start</th>
                <th>Dup End</th>
                <th>Size</th>
                <th>HR Status</th>
                <th>Microhomology</th>
                <th>Spanning Reads</th>
                <th>HQ Reads</th>
            </tr>
        </thead>
        <tbody>
            {''.join(rows)}
        </tbody>
    </table>
    <p class="table-note">Click any row for junction alignment detail.</p>
    """


def _build_module3_table(results):
    """Build summary table for module 3 results."""
    if not results:
        return "<p>No junctions quantified.</p>"

    rows = []
    for i, r in enumerate(results):
        jrr = r.get("JRR", 0)
        jrr_class = "jrr-high" if jrr > 0.1 else "jrr-med" if jrr > 0.01 else "jrr-low"

        row = f"""
        <tr class="clickable-row" data-target="detail-{i}">
            <td>{i + 1}</td>
            <td>{r.get('name', '')}</td>
            <td>{r.get('seq_id', '')}</td>
            <td>{r.get('dup_start', 0):,}</td>
            <td>{r.get('dup_end', 0):,}</td>
            <td>{r.get('dup_size', 0):,} bp</td>
            <td>{r.get('spanning_reads', 0)}</td>
            <td>{r.get('wt_mean_coverage', 0):.1f}x</td>
            <td class="{jrr_class}">{jrr:.4f}</td>
        </tr>"""
        rows.append(row)

    return f"""
    <div class="jrr-note">
        <strong>JRR (Junction Read Ratio)</strong> is proportional to duplication
        prevalence but is <em>not</em> an absolute allele frequency. It is best
        suited for comparing relative abundance across time-series samples.
    </div>
    <table class="results-table" id="summary-table">
        <thead>
            <tr>
                <th>#</th>
                <th>Name</th>
                <th>Sequence</th>
                <th>Dup Start</th>
                <th>Dup End</th>
                <th>Size</th>
                <th>Spanning Reads</th>
                <th>WT Coverage</th>
                <th>JRR</th>
            </tr>
        </thead>
        <tbody>
            {''.join(rows)}
        </tbody>
    </table>
    <p class="table-note">Click any row for junction detail.</p>
    """


# =============================================================================
# Detail panels (click-through views)
# =============================================================================

def _build_module1_details(results):
    """Build detail panels for module 1."""
    panels = []
    for i, r in enumerate(results):
        is_hr = r.get("is_hr_signature", False)
        hr_label = "HR" if is_hr else "non-HR"
        color = MECHANISM_COLORS.get("HR" if is_hr else "NHEJ", "#999")
        mh_seq = r.get("microhomology_seq", "")
        mh_display = f'<span class="mh-highlight">{mh_seq}</span>' if mh_seq else "none"

        hr_info = ""
        if is_hr:
            hr_info = f"""
                <div class="detail-item">
                    <label>HR Match Length</label>
                    <span>{r.get('hr_match_len', 0)} bp</span>
                </div>
                <div class="detail-item">
                    <label>HR Identity</label>
                    <span>{r.get('hr_identity', 0):.1%}</span>
                </div>
                <div class="detail-item">
                    <label>HR Scenario</label>
                    <span>{r.get('hr_scenario', '')}</span>
                </div>"""

        panel = f"""
        <div class="detail-panel" id="detail-{i}" style="display:none;">
            <h3>
                <span class="mechanism-tag" style="background:{color}">{hr_label}</span>
                Duplication #{i + 1}: {r.get('seq_id', '')}:{r.get('copy1_start', 0):,}-{r.get('copy2_end', 0):,}
            </h3>
            <div class="detail-grid">
                <div class="detail-item">
                    <label>Copy 1</label>
                    <span>{r.get('copy1_start', 0):,} - {r.get('copy1_end', 0):,}</span>
                </div>
                <div class="detail-item">
                    <label>Copy 2</label>
                    <span>{r.get('copy2_start', 0):,} - {r.get('copy2_end', 0):,}</span>
                </div>
                <div class="detail-item">
                    <label>Size</label>
                    <span>{r.get('size', 0):,} bp</span>
                </div>
                <div class="detail-item">
                    <label>Inter-copy Distance</label>
                    <span>{r.get('distance', 0):,} bp</span>
                </div>
                <div class="detail-item">
                    <label>Identity</label>
                    <span>{r.get('identity', 0):.1f}%</span>
                </div>
                <div class="detail-item">
                    <label>Orientation</label>
                    <span>{r.get('orientation', '')}</span>
                </div>
                <div class="detail-item">
                    <label>Microhomology</label>
                    <span>{r.get('microhomology_bp', 0)} bp: {mh_display}</span>
                </div>
                {hr_info}
            </div>
            <button class="back-btn" onclick="hideDetail({i})">← Back to summary</button>
        </div>"""
        panels.append(panel)

    return "\n".join(panels)


def _build_module2_details(results):
    """Build detail panels for module 2 with junction alignment."""
    panels = []
    for i, r in enumerate(results):
        mech = r.get("mechanism", "unknown")
        color = MECHANISM_COLORS.get(mech, "#999")
        junc_seq = r.get("junction_sequence", "")
        k = r.get("junction_k", 15)

        # Format junction sequence with midpoint marker
        if junc_seq:
            left = junc_seq[:k]
            right = junc_seq[k:]
            junc_display = (
                f'<span class="junc-left">{left}</span>'
                f'<span class="junc-break">|</span>'
                f'<span class="junc-right">{right}</span>'
            )
        else:
            junc_display = "N/A"

        panel = f"""
        <div class="detail-panel" id="detail-{i}" style="display:none;">
            <h3>
                <span class="mechanism-tag" style="background:{color}">{mech}</span>
                Junction #{i + 1}: {r.get('seq_id', '')}:{r.get('dup_start', 0):,}-{r.get('dup_end', 0):,}
            </h3>
            <div class="detail-grid">
                <div class="detail-item">
                    <label>Duplication</label>
                    <span>{r.get('dup_start', 0):,} - {r.get('dup_end', 0):,} ({r.get('dup_size', 0):,} bp)</span>
                </div>
                <div class="detail-item">
                    <label>Spanning Reads</label>
                    <span>{r.get('spanning_reads', 0)}</span>
                </div>
                <div class="detail-item">
                    <label>HQ Reads</label>
                    <span>{r.get('hq_reads', 0)}</span>
                </div>
                <div class="detail-item">
                    <label>Microhomology</label>
                    <span>{r.get('microhomology_bp', 0)} bp</span>
                </div>
            </div>
            <div class="junction-alignment">
                <h4>Junction Sequence (k={k})</h4>
                <div class="junction-legend">
                    <span class="junc-left-label">← End of copy 1</span>
                    <span class="junc-right-label">Start of copy 2 →</span>
                </div>
                <pre class="junction-seq">{junc_display}</pre>
            </div>
            <button class="back-btn" onclick="hideDetail({i})">← Back to summary</button>
        </div>"""
        panels.append(panel)

    return "\n".join(panels)


def _build_module3_details(results):
    """Build detail panels for module 3."""
    panels = []
    for i, r in enumerate(results):
        jrr = r.get("JRR", 0)

        panel = f"""
        <div class="detail-panel" id="detail-{i}" style="display:none;">
            <h3>Junction: {r.get('name', '')} ({r.get('seq_id', '')}:{r.get('dup_start', 0):,}-{r.get('dup_end', 0):,})</h3>
            <div class="detail-grid">
                <div class="detail-item">
                    <label>Duplication Size</label>
                    <span>{r.get('dup_size', 0):,} bp</span>
                </div>
                <div class="detail-item">
                    <label>Junction k</label>
                    <span>{r.get('junction_k', 0)} bp each side</span>
                </div>
                <div class="detail-item">
                    <label>Spanning Reads</label>
                    <span>{r.get('spanning_reads', 0)}</span>
                </div>
                <div class="detail-item">
                    <label>HQ Reads</label>
                    <span>{r.get('hq_reads', 0)}</span>
                </div>
                <div class="detail-item">
                    <label>WT Coverage</label>
                    <span>{r.get('wt_mean_coverage', 0):.1f}x</span>
                </div>
                <div class="detail-item jrr-detail">
                    <label>Junction Read Ratio (JRR)</label>
                    <span class="jrr-value">{jrr:.6f}</span>
                </div>
            </div>
            <div class="jrr-interpretation">
                <p>This JRR value indicates the relative abundance of this duplication
                junction in the population. Compare across time-series samples to
                track duplication dynamics.</p>
            </div>
            <button class="back-btn" onclick="hideDetail({i})">← Back to summary</button>
        </div>"""
        panels.append(panel)

    return "\n".join(panels)


# =============================================================================
# HTML page template
# =============================================================================

def _build_html_page(title, ref_name, table_html, detail_html,
                     module_name, extra_info):
    """Build complete HTML page."""
    module_labels = {
        "module1": "Reference Analysis",
        "module2": "Isolate Junction Discovery",
        "module3": "Population Quantification",
    }
    module_label = module_labels.get(module_name, module_name)

    extra_rows = ""
    if extra_info:
        for key, val in extra_info.items():
            extra_rows += f"<tr><td>{key}</td><td>{val}</td></tr>"

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            background: #f5f5f5;
            color: #333;
            line-height: 1.6;
        }}
        .header {{
            background: #2c3e50;
            color: white;
            padding: 20px 30px;
        }}
        .header h1 {{
            font-size: 24px;
            font-weight: 600;
        }}
        .header .subtitle {{
            font-size: 14px;
            opacity: 0.8;
            margin-top: 4px;
        }}
        .nav-bar {{
            background: #34495e;
            padding: 8px 30px;
            font-size: 13px;
            color: #bdc3c7;
        }}
        .nav-bar span {{
            margin-right: 20px;
        }}
        .nav-bar .active {{
            color: white;
            font-weight: 600;
        }}
        .container {{
            max-width: 1400px;
            margin: 20px auto;
            padding: 0 20px;
        }}
        .info-box {{
            background: white;
            border-radius: 6px;
            padding: 16px 20px;
            margin-bottom: 16px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }}
        .info-box table {{
            font-size: 13px;
        }}
        .info-box td {{
            padding: 3px 16px 3px 0;
        }}
        .info-box td:first-child {{
            font-weight: 600;
            color: #666;
        }}
        .results-table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            border-radius: 6px;
            overflow: hidden;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
            font-size: 13px;
        }}
        .results-table thead {{
            background: #ecf0f1;
        }}
        .results-table th {{
            padding: 10px 12px;
            text-align: left;
            font-weight: 600;
            color: #555;
            border-bottom: 2px solid #ddd;
        }}
        .results-table td {{
            padding: 8px 12px;
            border-bottom: 1px solid #eee;
        }}
        .results-table .clickable-row {{
            cursor: pointer;
        }}
        .results-table .clickable-row:hover {{
            background: #f0f7ff;
        }}
        .mechanism-tag {{
            display: inline-block;
            padding: 2px 8px;
            border-radius: 3px;
            color: white;
            font-size: 12px;
            font-weight: 600;
        }}
        .badge {{
            display: inline-block;
            padding: 1px 6px;
            border-radius: 3px;
            font-size: 11px;
            font-weight: 600;
        }}
        .badge-high {{ background: #d4edda; color: #155724; }}
        .badge-medium {{ background: #fff3cd; color: #856404; }}
        .badge-low {{ background: #f8d7da; color: #721c24; }}
        .table-note {{
            font-size: 12px;
            color: #888;
            margin-top: 8px;
            font-style: italic;
        }}
        .jrr-note {{
            background: #eef6ff;
            border-left: 4px solid #4A90D9;
            padding: 12px 16px;
            margin-bottom: 16px;
            font-size: 13px;
            border-radius: 0 4px 4px 0;
        }}
        .jrr-high {{ color: #c0392b; font-weight: 600; }}
        .jrr-med {{ color: #e67e22; font-weight: 600; }}
        .jrr-low {{ color: #27ae60; }}
        .detail-panel {{
            background: white;
            border-radius: 6px;
            padding: 24px;
            margin-top: 16px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }}
        .detail-panel h3 {{
            margin-bottom: 16px;
            font-size: 18px;
        }}
        .detail-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
            gap: 12px;
            margin-bottom: 20px;
        }}
        .detail-item {{
            background: #f8f9fa;
            padding: 10px 14px;
            border-radius: 4px;
        }}
        .detail-item label {{
            display: block;
            font-size: 11px;
            font-weight: 600;
            color: #888;
            text-transform: uppercase;
            margin-bottom: 4px;
        }}
        .detail-item span {{
            font-size: 14px;
        }}
        .junction-alignment {{
            background: #f8f9fa;
            border-radius: 4px;
            padding: 16px;
            margin-bottom: 16px;
        }}
        .junction-alignment h4 {{
            margin-bottom: 8px;
            font-size: 14px;
        }}
        .junction-legend {{
            display: flex;
            justify-content: space-between;
            font-size: 11px;
            color: #888;
            margin-bottom: 4px;
        }}
        .junction-seq {{
            font-family: "Courier New", monospace;
            font-size: 14px;
            letter-spacing: 1px;
            padding: 8px;
            background: white;
            border-radius: 3px;
            overflow-x: auto;
        }}
        .junc-left {{ color: #2980b9; }}
        .junc-right {{ color: #e74c3c; }}
        .junc-break {{ color: #333; font-weight: bold; margin: 0 2px; }}
        .mh-highlight {{
            background: #ffeaa7;
            padding: 1px 3px;
            border-radius: 2px;
            font-family: "Courier New", monospace;
        }}
        .jrr-detail {{ background: #eef6ff; }}
        .jrr-value {{ font-size: 20px; font-weight: 700; color: #2c3e50; }}
        .jrr-interpretation {{
            font-size: 13px;
            color: #666;
            padding: 12px;
            background: #f8f9fa;
            border-radius: 4px;
        }}
        .back-btn {{
            background: #ecf0f1;
            border: none;
            padding: 8px 16px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 13px;
            color: #555;
        }}
        .back-btn:hover {{ background: #ddd; }}
        .footer {{
            text-align: center;
            padding: 20px;
            font-size: 12px;
            color: #aaa;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Tandem</h1>
        <div class="subtitle">Detection and classification of tandem duplication mechanisms</div>
    </div>
    <div class="nav-bar">
        <span class="{'active' if module_name == 'module1' else ''}">Reference</span>
        <span class="{'active' if module_name == 'module2' else ''}">Isolate</span>
        <span class="{'active' if module_name == 'module3' else ''}">Population</span>
    </div>
    <div class="container">
        <div class="info-box">
            <table>
                <tr><td>Analysis Mode</td><td>{module_label}</td></tr>
                <tr><td>Reference</td><td>{ref_name}</td></tr>
                {extra_rows}
            </table>
        </div>

        <div id="summary-view">
            {table_html}
        </div>

        <div id="detail-views">
            {detail_html}
        </div>
    </div>

    <div class="footer">
        Generated by Tandem v0.2.0
    </div>

    <script>
        // Click-through from summary to detail
        document.querySelectorAll('.clickable-row').forEach(row => {{
            row.addEventListener('click', () => {{
                const target = row.getAttribute('data-target');
                document.getElementById('summary-view').style.display = 'none';
                document.getElementById(target).style.display = 'block';
            }});
        }});

        function hideDetail(idx) {{
            document.getElementById('detail-' + idx).style.display = 'none';
            document.getElementById('summary-view').style.display = 'block';
        }}

        // Table sorting
        document.querySelectorAll('.results-table th').forEach((th, colIdx) => {{
            th.style.cursor = 'pointer';
            th.addEventListener('click', () => {{
                const table = th.closest('table');
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));
                const isNumeric = rows.some(r => {{
                    const text = r.cells[colIdx]?.textContent.replace(/[,bp%x]/g, '').trim();
                    return !isNaN(parseFloat(text));
                }});
                const dir = th.dataset.sortDir === 'asc' ? 'desc' : 'asc';
                th.dataset.sortDir = dir;
                rows.sort((a, b) => {{
                    let aVal = a.cells[colIdx]?.textContent.replace(/[,bp%x]/g, '').trim() || '';
                    let bVal = b.cells[colIdx]?.textContent.replace(/[,bp%x]/g, '').trim() || '';
                    if (isNumeric) {{
                        aVal = parseFloat(aVal) || 0;
                        bVal = parseFloat(bVal) || 0;
                        return dir === 'asc' ? aVal - bVal : bVal - aVal;
                    }}
                    return dir === 'asc' ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
                }});
                rows.forEach(r => tbody.appendChild(r));
            }});
        }});
    </script>
</body>
</html>"""
