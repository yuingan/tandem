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
    ok = len(found) == len(multiples) >= 1 and all(c["unit"] == 113 for c in multiples)
    results.append(("multiples of a 113 bp unit, --min-size 200", ok,
                    f"{len(multiples)} of {len(found)} set aside as multiples"))
    found = calls(rows, min_size=100)
    kept = [c for c in found if c["structure"] == "tandem"]
    ok = (len(kept) == 1 and kept[0]["size"] == 113 and kept[0]["copy_number"] == 6
          and all(c["structure"] == "unit_multiple" for c in found if c not in kept))
    detail = (f"{kept[0]['size']} bp x{kept[0]['copy_number']}"
              if len(kept) == 1 else f"{len(kept)} tandem call(s)")
    results.append(("the same at --min-size 100", ok, detail))

    # P. aeruginosa 2,036,966-2,042,137, the show-coords rows from the dev8
    # run: about 21 copies of a 246 bp unit at 80-88% identity. The
    # alignment at 246 bp breaks into two pieces; those at 492, 738 and 984
    # bp overlap the first; pairs at 984-2,214 bp span both. One array, and
    # every pair in it a multiple.
    rows = [(2036966, 2037445, 2037224, 2037703, 85.42), (2037731, 2038871, 2037977, 2039117, 82.29),
            (2037933, 2039117, 2038425, 2039609, 82.30), (2037934, 2039638, 2038918, 2040622, 81.06),
            (2037960, 2039857, 2038698, 2040595, 83.73), (2038688, 2039872, 2040902, 2042086, 80.25),
            (2038698, 2039609, 2039928, 2040839, 80.24), (2038960, 2039869, 2040928, 2041837, 80.59),
            (2039717, 2040625, 2040947, 2041855, 81.31), (2039921, 2040625, 2040905, 2041609, 82.68),
            (2039930, 2040415, 2041652, 2042137, 83.57), (2040161, 2040376, 2040899, 2041114, 88.53),
            (2040885, 2041891, 2041131, 2042137, 84.14)]
    found = calls([aln(*r) for r in rows] + [aln(r[2], r[3], r[0], r[1], r[4]) for r in rows])
    arrays = [c for c in found if c["structure"] == "tandem"]
    left = [c for c in found if c["structure"] not in ("tandem", "unit_multiple")]
    ok = (len(arrays) == 1 and abs(arrays[0]["size"] - 246) <= 12
          and 20 <= arrays[0]["copy_number"] <= 22 and not left)
    detail = (f"{arrays[0]['size']} bp x{arrays[0]['copy_number']}, {len(left)} other call(s)"
              if len(arrays) == 1 else f"{len(arrays)} tandem call(s), {len(left)} other")
    results.append(("degenerate array, P. aeruginosa 2.04 Mb", ok, detail))

    # Shaped like E. coli 1,466,052: 3.2 copies of a 294 bp unit, so the
    # alignment at 294 bp covers 2.2 units, and the pair at 588 bp is its
    # multiple. Arrays end in partial units; this one used to vanish.
    found = calls([aln(1466052, 1466711, 1466346, 1467005, 89.5),
                   aln(1466052, 1466417, 1466640, 1467005, 89.1)])
    arrays = [c for c in found if c["structure"] == "tandem"]
    ok = (len(arrays) == 1 and arrays[0]["size"] == 294
          and arrays[0]["copy_number"] == 3.2 and len(found) == 2)
    detail = (f"{arrays[0]['size']} bp x{arrays[0]['copy_number']}"
              if len(arrays) == 1 else f"{len(arrays)} tandem call(s)")
    results.append(("array ending in a partial unit", ok, detail))

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
