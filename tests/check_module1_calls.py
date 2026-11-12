"""How Module 1 turns self-alignments into calls.

    python tests/check_module1_calls.py

Each case rebuilds, from the coordinates in the species runs, the alignments
behind a call that came out wrong, feeds them to identify_tandem_duplications
and checks what comes out.
"""

import sys

from tandem.module1_reference import identify_tandem_duplications


def aln(ref_start, ref_end, qry_start, qry_end, identity=95.0, seq="chr"):
    return {"ref_id": seq, "qry_id": seq, "identity": identity,
            "ref_start": ref_start, "ref_end": ref_end,
            "qry_start": qry_start, "qry_end": qry_end,
            "ref_alen": abs(ref_end - ref_start) + 1,
            "qry_alen": abs(qry_end - qry_start) + 1}


def calls(alignments, length=7_000_000, **kw):
    return identify_tandem_duplications(alignments, {"chr": length}, **kw)


def main():
    results = []

    # SBW25 6,485,658: one pair at 88% identity, aligned from each side and
    # stopping at different bases. Was reported twice.
    found = calls([aln(6485658, 6485917, 6485976, 6486229, 88.5),
                   aln(6485976, 6486283, 6485658, 6485971, 88.9)])
    ok = len(found) == 1 and found[0]["copy2_end"] == 6486283
    results.append(("one pair aligned from both sides", ok,
                    f"{len(found)} call(s)"))

    # Shaped like E. coli 2,304,390, where an array at 226 bp and a pair at
    # 339 bp point to six copies of a 113 bp unit.
    rows = [aln(2304390, 2304954, 2304503, 2305067),
            aln(2304390, 2304841, 2304616, 2305067),
            aln(2304390, 2304728, 2304729, 2305067)]
    found = calls(rows)
    multiples = [c for c in found if c["structure"] == "unit_multiple"]
    ok = len(found) == 2 and len(multiples) == 2 and all(c["unit"] == 113 for c in multiples)
    results.append(("multiples of a 113 bp unit, --min-size 200", ok,
                    f"{len(multiples)} of {len(found)} set aside as multiples"))
    found = calls(rows, min_size=100)
    kept = [c for c in found if c["structure"] == "tandem"]
    ok = (len(found) == 3 and len(kept) == 1
          and kept[0]["size"] == 113 and kept[0]["copy_number"] == 6)
    detail = (f"{kept[0]['size']} bp x{kept[0]['copy_number']}"
              if len(kept) == 1 else f"{len(kept)} tandem call(s)")
    results.append(("the same at --min-size 100", ok, detail))

    # Alignments whose shift is not the same at both ends: one starts on the
    # diagonal, the other drifts 2 bp to 300 bp through a low-complexity
    # stretch. Neither is a repeat unit. dev7 divided by the first one's zero
    # shift, and read the second as a 2 bp unit and dropped the pair it
    # covers. The pair is E. coli 1,096,961 (period 362), one of three calls
    # dev7 lost there; the drifting alignments are made up.
    found = calls([aln(1096961, 1097221, 1097323, 1097583, 99.2),
                   aln(1096900, 1097700, 1096900, 1098000, 90.0),
                   aln(1096900, 1097700, 1096902, 1098000, 90.0)])
    ok = len(found) == 1 and found[0]["structure"] == "dispersed_pair"
    results.append(("alignments drifting off the diagonal", ok,
                    found[0]["structure"] if len(found) == 1 else f"{len(found)} call(s)"))

    # S. coelicolor: the terminal inverted repeats meet through the origin in
    # circular mode, one base apart. Was reported as a 21 kb tandem.
    found = calls([aln(1, 21653, 8667507, 8645855, 100.0)], length=8667507)
    ok = len(found) == 1 and found[0]["structure"] == "inverted_pair"
    results.append(("terminal inverted repeats, circular", ok,
                    found[0]["structure"] if found else "no call"))

    # Four copies of a 1.2 kb element, evenly spaced. Pairs at the same
    # spacing are still different pairs.
    starts = [5000, 12200, 19400, 26600]
    found = calls([aln(a, a + 1199, b, b + 1199, 99.0)
                   for a in starts for b in starts if a < b])
    results.append(("evenly spaced family of four", len(found) == 6,
                    f"{len(found)} pair(s)"))

    # A plain two-copy tandem is untouched.
    found = calls([aln(20000, 24999, 25000, 29999, 100.0),
                   aln(25000, 29999, 20000, 24999, 100.0)])
    ok = len(found) == 1 and found[0]["structure"] == "tandem"
    results.append(("two-copy tandem", ok, f"{len(found)} call(s)"))

    failed = 0
    for name, ok, detail in results:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<42} {detail}")
        failed += not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
