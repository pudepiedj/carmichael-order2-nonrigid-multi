#!/usr/bin/env python3
"""cn2_multi_afile.py -- write the OEIS a-file data/a175531_multi.txt from data/cn2_nonrigid_multi.tsv.

    python3 tools/cn2_multi_afile.py            (run from the repository root)

Part 1 lists the numbers with three primes p at which k == p (mod p^2-1), Part 2 those with two; each part is sorted.
"""
import csv

SRC, OUT = 'data/cn2_nonrigid_multi.tsv', 'data/a175531_multi.txt'

rows = list(csv.DictReader(open(SRC), delimiter='\t'))
parts = {k: [r for r in rows if r['type_p'].count(',') == k - 1] for k in (3, 2)}
for k in parts:
    parts[k].sort(key=lambda r: int(r['n']))

head = f"""\
# Terms of A175531 with k == p (mod p^2-1) at two or three primes p | k (and k == 1 (mod q^2-1) at the others).
# These are EXAMPLES, not a complete list: they were found by Howe's construction on particular moduli, not by
# searching all k below a bound. Every number below was verified directly: k is squarefree, and for every prime
# p | k, k == 1 or k == p (mod p^2-1). The third column lists the p with k == p (mod p^2-1), joined by '.'.
# {len(parts[2])} numbers with two such primes and {len(parts[3])} with three. Earlier examples of both kinds are in the
# a-file of non-rigid examples already linked from A175531; all of those are included here.
# Source and verifier: https://github.com/pudepiedj/carmichael-order2-nonrigid-multi
#"""
with open(OUT, 'w') as o:
    o.write(head + '\n')
    for i, (k, label) in enumerate(((3, 'three'), (2, 'two')), 1):
        o.write(f"# Part {i}: the {len(parts[k])} terms found with {label} such primes p.\n")
        o.write("# Columns: k  log10(k)  primes p with k == p (mod p^2-1)  prime factors\n")
        for r in parts[k]:
            o.write(f"{r['n']} {float(r['log10']):.2f} {r['type_p'].replace(',', '.')} {r['primes']}\n")
        o.write("#\n")
print(f"wrote {OUT}: {len(parts[3])} + {len(parts[2])} terms")
