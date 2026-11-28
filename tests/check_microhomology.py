"""Microhomology at a junction, on sequences where the answer is planted.

    python tests/check_microhomology.py

The junction joins the end of the unit to its start. Microhomology is the
bases there that could have come from either side: the end of the unit read
against the bases just before it, and its start against the bases just after
the second copy. Checks both measures in tandem/junction.py, and the case
that used to go wrong: a unit that begins and ends with the same bases.
"""

import random
import sys

from tandem.junction import junction_microhomology, single_copy_microhomology


def rnd(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def main():
    random.seed(4)
    results = []

    # X ends with ACG like the unit does (3 bp to the left); Y starts with TT
    # like the unit does (2 bp to the right). 5 bp in all.
    X = rnd(497) + "ACG"
    Y = "TT" + rnd(498)
    U = "TT" + rnd(296) + "ACG"
    while X[-4] == U[-4] or Y[2] == U[2]:
        U = "TT" + rnd(296) + "ACG"

    s0, e0 = len(X), len(X) + len(U)
    n, seq = single_copy_microhomology(X + U + Y, s0, e0)
    results.append(("one reference copy", n == 5 and seq == "ACGTT", f"{n} bp {seq}"))

    product = X + U + U + Y
    a, p = len(X), len(U)
    n, seq = junction_microhomology(product, a, a + p, a + p, a + 2 * p)
    results.append(("two adjacent copies", n == 5 and seq == "ACGTT", f"{n} bp {seq}"))

    # as NUCmer reports it: the copies overlap by the microhomology, then are
    # made adjacent with the first copy keeping the overlap
    n, seq = junction_microhomology(product, a - 3, a + p + 2, a + p + 2,
                                    a + 2 * p + 2, overlap=5)
    results.append(("copies after overlap resolution", n == 5, f"{n} bp"))

    # Unit begins and ends with GTG, but neither flank matches: no
    # microhomology. The old measure called this 3 bp.
    X = rnd(499) + "T"
    Y = "C" + rnd(499)
    U = "GTGA" + rnd(292) + "GGTG"
    n, _ = single_copy_microhomology(X + U + Y, len(X), len(X) + len(U))
    results.append(("unit starting and ending alike", n == 0, f"{n} bp"))

    failed = 0
    for name, ok, detail in results:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<42} {detail}")
        failed += not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
