#!/usr/bin/env python3
"""Generate a random PFSP instance in the layout of testing/DD_Ta055.txt.

usage: gen_instance.py <njobs> <nmac> <seed> <out>

Layout (see PfspInstance::readDataFromFile): header "njobs nmac"; one line per
job of "machine_index time" pairs (times uniform 1..99); a "Reldue" keyword
line; then one line per job "-1 duedate -1 priority".  Due dates are uniform
in [sum of the job's times, 4 * that sum], priorities uniform 1..10.
"""
import random
import sys


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    n, m, seed, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    rng = random.Random(seed)
    times = [[rng.randint(1, 99) for _ in range(m)] for _ in range(n)]
    with open(out, "w") as f:
        f.write("%d %d\n" % (n, m))
        for row in times:
            f.write(" ".join("%3d %2d" % (k, t) for k, t in enumerate(row)) + "\n")
        f.write("Reldue\n")
        for row in times:
            s = sum(row)
            f.write("-1\t%d\t-1\t%d\n" % (rng.randint(s, 4 * s), rng.randint(1, 10)))


if __name__ == "__main__":
    main()
