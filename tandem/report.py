"""HTML report for Tandem results.

A summary table of the duplications or junctions found, each row opening a
detail view, in the spirit of breseq's output.
"""

import logging
import os
from pathlib import Path

from . import __version__

logger = logging.getLogger("tandem")


ARCHITECTURE_COLORS = {
    "flanking_repeat": "#7ED321",
    "no_flanking_repeat": "#4A90D9",
    "array": "#F5A623",
    "dispersed_pair": "#999999",
}

# results of the legacy search still carry a mechanism label
MECHANISM_COLORS = {
    "NHEJ": "#4A90D9",
    "MMEJ": "#F5A623",
    "SSA": "#D0021B",
    "HR": "#7ED321",
    "unknown": "#999999",
}

SHOW_SEQUENCES_UP_TO = 2000    # bp per copy; longer copies are only in copies.fasta
SHOW_COPIES_UP_TO = 30


# =============================================================================
# Report generation
# =============================================================================

def generate_report(results, module_name, output_dir, ref_name="",
                    extra_info=None):
    """Generate an HTML report for Tandem results."""
    output_dir = Path(output_dir)
    report_path = output_dir / "tandem_report.html"

    if module_name == "module1":
        title = "Tandem - Reference Genome Analysis"
        table_html = _build_module1_table(results)
        detail_html = _build_module1_details(results, _read_fasta(output_dir / "copies.fasta"))
    elif module_name == "module2":
        title = "Tandem - Isolate Junction Discovery"
        if results and "spanning_reads" in results[0]:
            table_html = _build_legacy_module2_table(results)
            detail_html = _build_legacy_module2_details(results)
        else:
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


def _read_fasta(path):
    if not os.path.isfile(path):
        return {}
    seqs, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def _num(value, digits=4, percent=False):
    if value is None or value == "":
        return "–"
    return f"{value:.0%}" if percent else f"{value:.{digits}f}"


def _item(label, value):
    return f"""
                <div class="detail-item">
                    <label>{label}</label>
                    <span>{value}</span>
                </div>"""


def _architecture(r):
    if r.get("architecture"):
        return r["architecture"]
    return "flanking_repeat" if r.get("flanking_repeat") else "no_flanking_repeat"


def _tag(architecture):
    color = ARCHITECTURE_COLORS.get(architecture, "#999")
    return (f'<span class="mechanism-tag" style="background:{color}">'
            f'{architecture.replace("_", " ")}</span>')


def _repeat_position(geometry):
    """Where the flanking repeat sits, in words, from junction.py's labels."""
    g = geometry or ""
    if g.startswith("s1"):
        return "at the start of each copy and after copy 2"
    if g.startswith("s2"):
        return "at the end of each copy and before copy 1"
    for code, words in (("_s3", "outside both ends"),
                        ("_s1", "outside the start, inside the end"),
                        ("_s2", "inside the start, outside the end"),
                        ("_s4", "inside both ends")):
        if code in g:
            return words
    return g.replace("_", " ")


def _junction_display(junction_sequence, mh=0):
    """Unit end | unit start, with the microhomology marked on both sides."""
    if not junction_sequence or "|" not in junction_sequence:
        return "N/A"
    left, right = junction_sequence.split("|", 1)
    if mh:
        right = f'<span class="mh-highlight">{right[:mh]}</span>{right[mh:]}'
    return (f'<span class="junc-left">{left}</span>'
            f'<span class="junc-break">|</span>'
            f'<span class="junc-right">{right}</span>')


# =============================================================================
# Module 1
# =============================================================================

def _build_module1_table(results):
    """Build summary table for module 1 results."""
    if not results:
        return "<p>No tandem duplications detected.</p>"

    rows = []
    for i, r in enumerate(results):
        architecture = _architecture(r)
        mh = "–" if architecture == "array" else f"{r.get('microhomology_bp', 0)} bp"
        rows.append(f"""
        <tr class="clickable-row" data-target="detail-{i}">
            <td>{i + 1}</td>
            <td>{r.get('seq_id', '')}</td>
            <td>{r.get('copy1_start', 0):,}</td>
            <td>{r.get('copy1_end', 0):,}</td>
            <td>{r.get('copy2_start', 0):,}</td>
            <td>{r.get('copy2_end', 0):,}</td>
            <td>{r.get('size', 0):,} bp</td>
            <td>{r.get('copy_number', 2)}</td>
            <td>{r.get('distance', 0):,} bp</td>
            <td>{r.get('identity', 0):.1f}%</td>
            <td>{_tag(architecture)}</td>
            <td>{mh}</td>
        </tr>""")

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
                <th>Copy Size</th>
                <th>Copies</th>
                <th>Gap</th>
                <th>Identity</th>
                <th>Junction</th>
                <th>Microhomology</th>
            </tr>
        </thead>
        <tbody>
            {''.join(rows)}
        </tbody>
    </table>
    <p class="table-note">For arrays, copies 1 and 2 are the first two units. Click
    any row for details and the copy sequences.</p>
    """


NOTES = {
    "copies_diverged": "Copies below the --diverged-below identity: the junction "
                       "shows the sequence as it is now, not as the duplication left it.",
    "junction_distant": "Copies far apart (beyond --flag-distance).",
    "too_small_for_repeat_search": "Copies too short to search for a flanking repeat.",
}


def _build_module1_details(results, copies):
    """Build detail panels for module 1."""
    panels = []
    for i, r in enumerate(results):
        architecture = _architecture(r)
        n = r.get("copy_number", 2)
        unit = r.get("size", 0)

        if architecture == "array":
            end = r["copy1_start"] + round(n * unit) - 1
            span = f"{r.get('seq_id', '')}:{r['copy1_start']:,}-{end:,}"
            items = [
                _item("Unit", f"{unit:,} bp"),
                _item("Copies", n),
                _item("Array", f"{r['copy1_start']:,} - {end:,}"),
                _item("Identity one unit apart", f"{r.get('identity', 0):.1f}%"),
            ]
            junction_note = ("<p>Each junction inside an array has a full copy of the "
                             "unit on either side, so its architecture says nothing "
                             "about how the array formed.</p>")
        else:
            span = f"{r.get('seq_id', '')}:{r['copy1_start']:,}-{r['copy2_end']:,}"
            mh_seq = r.get("microhomology_seq", "")
            mh = f'<span class="mh-highlight">{mh_seq}</span>' if mh_seq else "none"
            items = [
                _item("Copy 1", f"{r['copy1_start']:,} - {r['copy1_end']:,}"),
                _item("Copy 2", f"{r['copy2_start']:,} - {r['copy2_end']:,}"),
                _item("Copy size", f"{unit:,} bp"),
                _item("Gap between copies", f"{r.get('distance', 0):,} bp"),
                _item("Identity", f"{r.get('identity', 0):.1f}%"),
                _item("Orientation", r.get("orientation", "")),
                _item("Microhomology", f"{r.get('microhomology_bp', 0)} bp: {mh}"),
            ]
            if r.get("flanking_repeat"):
                items += [
                    _item("Flanking repeat", f"{r.get('repeat_len', 0)} bp"),
                    _item("Repeat identity", f"{r.get('repeat_identity', 0):.1%}"),
                    _item("Repeat position", _repeat_position(r.get("repeat_geometry"))),
                ]
            junction_note = ""

        notes = [NOTES.get(x, x) for x in str(r.get("classification_note", "")).split(";") if x]
        note_html = "".join("<p>" + x + "</p>" for x in notes)
        items_html = "".join(items)
        tag = _tag(architecture)
        sequences = _copy_sequences(i + 1, copies, unit, architecture == "array")

        panels.append(f"""
        <div class="detail-panel" id="detail-{i}" style="display:none;">
            <h3>
                {tag}
                Duplication #{i + 1}: {span}
            </h3>
            <div class="detail-grid">{items_html}
            </div>
            <div class="detail-notes">{junction_note}{note_html}</div>
            {sequences}
            <button class="back-btn" onclick="hideDetail({i})">← Back to summary</button>
        </div>""")

    return "\n".join(panels)


def _copy_sequences(row, copies, unit, is_array):
    """Each copy's sequence, bases differing from copy 1 marked."""
    names = [k for k in copies if k.split(".")[0] == str(row)]
    if not names:
        return ""
    if unit > SHOW_SEQUENCES_UP_TO or len(names) > SHOW_COPIES_UP_TO:
        return (f'<div class="copy-block"><h4>Copy sequences</h4><p>{len(names)} copies '
                f'of {unit:,} bp: records {row}.1 to {row}.{len(names)} in '
                f'copies.fasta.</p></div>')

    first = copies[names[0]]
    blocks = []
    for name in names:
        seq = copies[name]
        # base by base, no gaps: only meaningful for copies of the same length
        # (or an array's last, partial unit)
        if name != names[0] and (len(seq) == len(first) or (is_array and len(seq) < len(first))):
            seq = "".join(b if b == a else '<span class="diff">' + b + "</span>"
                          for a, b in zip(first, seq))
        number = name.split(".")[1]
        blocks.append(f'<div class="copy-name">Copy {number} ({name} in copies.fasta)</div>'
                      f'<div class="copy-seq">{seq}</div>')
    return ('<div class="copy-block"><h4>Copy sequences</h4>'
            '<p class="table-note">Bases that differ from copy 1 are marked, position '
            'by position; insertions and deletions are not aligned.</p>'
            + "".join(blocks) + "</div>")


# =============================================================================
# Module 2
# =============================================================================

def _build_module2_table(results):
    """Junctions found by the scan."""
    if not results:
        return "<p>No junctions found.</p>"

    given = "provided_start" in results[0]
    control = "control_reads" in results[0]
    rows = []
    for i, r in enumerate(results):
        cells = [
            f"<td>{i + 1}</td>",
            f"<td>{r.get('seq_id', '')}</td>",
        ]
        if given:
            cells.append(f"<td>{r['provided_start']:,} - {r['provided_end']:,}</td>")
        cells += [
            f"<td>{r['start']:,}</td>",
            f"<td>{r['end']:,}</td>",
            f"<td>{r['size']:,} bp</td>",
            f"<td>{r.get('microhomology', 0)} bp</td>",
            f"<td>{r['junction_reads']}</td>",
            f"<td>{r['reference_reads']}</td>",
            f"<td>{_num(r.get('dosage'), 3)}</td>",
            f"<td>{r.get('dosage_from', '')}</td>",
            f"<td>{_num(r.get('depth_dosage'), 3)}</td>",
        ]
        if control:
            cells.append(f"<td>{r.get('control_reads', 0)}</td>")
        cells_html = "".join(cells)
        rows.append(f"""
        <tr class="clickable-row" data-target="detail-{i}">
            {cells_html}
        </tr>""")

    headers = ["#", "Sequence"] + (["Given"] if given else []) + [
        "Start", "End", "Size", "Microhomology", "Junction Reads",
        "Reference Reads", "Dosage", "From", "Depth Estimate",
    ] + (["Control Reads"] if control else [])
    header_html = "".join("<th>" + h + "</th>" for h in headers)

    return f"""
    <div class="jrr-note">
        Start and end are the first and last base of the duplicated unit.
        <strong>Dosage</strong> is junction reads over reference reads, each
        corrected for the read placements that can be counted there: the extra
        copies per cell, 1 for a duplication every cell carries. Where no read
        can be placed across the junction, it comes from depth instead.
    </div>
    <table class="results-table" id="summary-table">
        <thead>
            <tr>{header_html}</tr>
        </thead>
        <tbody>
            {''.join(rows)}
        </tbody>
    </table>
    <p class="table-note">Click any row for the junction sequence.</p>
    """


def _build_module2_details(results):
    panels = []
    for i, r in enumerate(results):
        items = [
            _item("Duplicated unit", f"{r['start']:,} - {r['end']:,} ({r['size']:,} bp)"),
            _item("Junction reads", r["junction_reads"]),
            _item("Reference reads", r["reference_reads"]),
            _item("Dosage", f"{_num(r.get('dosage'))} (from {r.get('dosage_from') or '–'})"),
            _item("Dosage before correction", _num(r.get("dosage_raw"))),
            _item("Dosage from depth", _num(r.get("depth_dosage"))),
            _item("Junction countable", _num(r.get("junction_visible"), percent=True)),
            _item("Reference countable", _num(r.get("reference_visible"), percent=True)),
        ]
        if "provided_start" in r:
            items.insert(1, _item("Given", f"{r['provided_start']:,} - {r['provided_end']:,}"))
        if "control_reads" in r:
            items.append(_item("Reads in the control", r["control_reads"]))
        items_html = "".join(items)
        junction = _junction_display(r.get("junction_sequence"), r.get("microhomology", 0))

        panels.append(f"""
        <div class="detail-panel" id="detail-{i}" style="display:none;">
            <h3>Junction #{i + 1}: {r.get('seq_id', '')}:{r['start']:,}-{r['end']:,}</h3>
            <div class="detail-grid">{items_html}
            </div>
            <div class="junction-alignment">
                <h4>Junction sequence</h4>
                <div class="junction-legend">
                    <span>← end of the unit (copy 1)</span>
                    <span>start of the unit (copy 2) →</span>
                </div>
                <pre class="junction-seq">{junction}</pre>
                <p class="table-note">Marked: microhomology, the bases the junction
                shares with the reference on both sides; the junction can be placed
                anywhere along it.</p>
            </div>
            <div class="detail-notes">
                <p>Countable: the share of read placements across the junction, or
                across the reference boundary, with a unique 20-mer on both sides.
                Below 100% a repeat lies near the junction and the dosage is
                corrected for it.</p>
            </div>
            <button class="back-btn" onclick="hideDetail({i})">← Back to summary</button>
        </div>""")

    return "\n".join(panels)


def _build_legacy_module2_table(results):
    """Results of --legacy-search."""
    rows = []
    for i, r in enumerate(results):
        is_hr = r.get("flanking_repeat", False)
        architecture = "flanking_repeat" if is_hr else "no_flanking_repeat"
        rows.append(f"""
        <tr class="clickable-row" data-target="detail-{i}">
            <td>{i + 1}</td>
            <td>{r.get('seq_id', '')}</td>
            <td>{r.get('dup_start', 0):,}</td>
            <td>{r.get('dup_end', 0):,}</td>
            <td>{r.get('dup_size', 0):,} bp</td>
            <td>{_tag(architecture)}</td>
            <td>{r.get('microhomology_bp', 0)} bp</td>
            <td>{r.get('spanning_reads', 0)}</td>
            <td>{r.get('hq_reads', 0)}</td>
        </tr>""")

    return f"""
    <table class="results-table" id="summary-table">
        <thead>
            <tr>
                <th>#</th>
                <th>Sequence</th>
                <th>Dup Start</th>
                <th>Dup End</th>
                <th>Size</th>
                <th>Flanking repeat</th>
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


def _build_legacy_module2_details(results):
    panels = []
    for i, r in enumerate(results):
        mech = r.get("mechanism", "unknown")
        color = MECHANISM_COLORS.get(mech, "#999")
        junc_seq = r.get("junction_sequence", "")
        k = r.get("junction_k", 15)
        display = _junction_display(junc_seq[:k] + "|" + junc_seq[k:]) if junc_seq else "N/A"
        items_html = "".join([
            _item("Duplication", f"{r.get('dup_start', 0):,} - {r.get('dup_end', 0):,} "
                                 f"({r.get('dup_size', 0):,} bp)"),
            _item("Spanning reads", r.get("spanning_reads", 0)),
            _item("HQ reads", r.get("hq_reads", 0)),
            _item("Microhomology", f"{r.get('microhomology_bp', 0)} bp"),
        ])

        panels.append(f"""
        <div class="detail-panel" id="detail-{i}" style="display:none;">
            <h3>
                <span class="mechanism-tag" style="background:{color}">{mech}</span>
                Junction #{i + 1}: {r.get('seq_id', '')}:{r.get('dup_start', 0):,}-{r.get('dup_end', 0):,}
            </h3>
            <div class="detail-grid">{items_html}
            </div>
            <div class="junction-alignment">
                <h4>Junction Sequence (k={k})</h4>
                <div class="junction-legend">
                    <span>← End of copy 1</span>
                    <span>Start of copy 2 →</span>
                </div>
                <pre class="junction-seq">{display}</pre>
            </div>
            <button class="back-btn" onclick="hideDetail({i})">← Back to summary</button>
        </div>""")

    return "\n".join(panels)


# =============================================================================
# Module 3
# =============================================================================

def _build_module3_table(results):
    """Build summary table for module 3 results."""
    if not results:
        return "<p>No junctions quantified.</p>"

    rows = []
    for i, r in enumerate(results):
        dosage = r.get("dosage")
        level = dosage or 0
        dosage_class = "jrr-high" if level > 0.1 else "jrr-med" if level > 0.01 else "jrr-low"

        rows.append(f"""
        <tr class="clickable-row" data-target="detail-{i}">
            <td>{i + 1}</td>
            <td>{r.get('name', '')}</td>
            <td>{r.get('seq_id', '')}</td>
            <td>{r.get('dup_start', 0):,}</td>
            <td>{r.get('dup_end', 0):,}</td>
            <td>{r.get('dup_size', 0):,} bp</td>
            <td>{r.get('junction_reads', 0)}</td>
            <td>{r.get('reference_reads', 0)}</td>
            <td class="{dosage_class}">{_num(dosage)}</td>
            <td>{r.get('dosage_from', '')}</td>
            <td>{_num(r.get('depth_dosage'))}</td>
        </tr>""")

    return f"""
    <div class="jrr-note">
        <strong>Dosage</strong> is junction reads over reference reads at the
        same boundary, each corrected for the read placements that can be
        counted there. It estimates the extra copies per cell; for a duplication
        adding one copy, the share of cells that carry it. Where no read can be
        placed across the junction, it comes from depth instead.
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
                <th>Junction Reads</th>
                <th>Reference Reads</th>
                <th>Dosage</th>
                <th>From</th>
                <th>Depth Estimate</th>
            </tr>
        </thead>
        <tbody>
            {''.join(rows)}
        </tbody>
    </table>
    <p class="table-note">Click any row for junction detail.</p>
    """


def _build_module3_details(results):
    """Build detail panels for module 3."""
    panels = []
    for i, r in enumerate(results):
        moved = ""
        if (r.get("dup_start"), r.get("dup_end")) != (r.get("provided_start"), r.get("provided_end")):
            moved = (f"<p>Given as {r.get('provided_start', 0):,} - {r.get('provided_end', 0):,}; "
                     f"measured at the junction the reads show.</p>")
        items = [
            _item("Duplication size", f"{r.get('dup_size', 0):,} bp"),
            _item("Junction reads", r.get("junction_reads", 0)),
            _item("Reference reads", r.get("reference_reads", 0)),
            _item("Dosage before correction", _num(r.get("dosage_raw"))),
            _item("Dosage from depth", _num(r.get("depth_dosage"))),
            _item("Junction countable", _num(r.get("junction_visible"), percent=True)),
            _item("Microhomology", f"{r.get('microhomology_bp', 0)} bp"),
            _item("Flanking repeat", f"{r.get('repeat_len', 0)} bp" if r.get("flanking_repeat") else "none"),
        ]
        items_html = "".join(items)
        source = r.get("dosage_from") or "–"
        dosage = _num(r.get("dosage"))
        panels.append(f"""
        <div class="detail-panel" id="detail-{i}" style="display:none;">
            <h3>Junction: {r.get('name', '')} ({r.get('seq_id', '')}:{r.get('dup_start', 0):,}-{r.get('dup_end', 0):,})</h3>
            <div class="detail-grid">{items_html}
                <div class="detail-item jrr-detail">
                    <label>Dosage (from {source})</label>
                    <span class="jrr-value">{dosage}</span>
                </div>
            </div>
            <div class="jrr-interpretation">
                {moved}
                <p>Extra copies per cell. For a duplication adding one copy, this
                is the share of cells that carry it.</p>
            </div>
            <button class="back-btn" onclick="hideDetail({i})">← Back to summary</button>
        </div>""")

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
            overflow-wrap: anywhere;
        }}
        .detail-notes {{
            font-size: 13px;
            color: #666;
            margin-bottom: 16px;
        }}
        .copy-block {{
            background: #f8f9fa;
            border-radius: 4px;
            padding: 16px;
            margin-bottom: 16px;
        }}
        .copy-block h4 {{
            margin-bottom: 8px;
            font-size: 14px;
        }}
        .copy-block p {{
            font-size: 13px;
            color: #666;
        }}
        .copy-name {{
            font-size: 12px;
            color: #666;
            margin-top: 10px;
        }}
        .copy-seq {{
            font-family: "Courier New", monospace;
            font-size: 12px;
            background: white;
            padding: 6px 8px;
            border-radius: 3px;
            word-break: break-all;
            max-height: 14em;
            overflow-y: auto;
        }}
        .diff {{ background: #fab1a0; font-weight: 700; }}
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
            word-break: break-all;
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
        <div class="subtitle">Tandem duplication detection and junction architecture</div>
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
        Generated by Tandem v{__version__}
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
