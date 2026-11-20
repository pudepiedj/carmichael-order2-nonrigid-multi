#!/usr/bin/env python3
"""Build the composite OEIS a-file for A175531 (examples, not complete lists).

Reads the three datasets of carmichael-order2-nonrigid-multi (data/cn2_nonrigid_multi.tsv,
data/cn2_ktree_multi.tsv, data/cn2_divisible_5_7.tsv) and writes one a-file in five parts:
  1. terms divisible by 5
  2. terms divisible by 7
  3. four type-p primes
  4. three type-p primes
  5. two type-p primes

A type-p prime of k is a prime p | k with k == p (mod p^2-1); every other prime q | k has
k == 1 (mod q^2-1).  After writing, the output file is re-read and every row is checked from
scratch: distinct primes, product = k, every prime satisfies one of the two congruences, and the
listed type-p set is exactly the set with k == p.

    python3 cn2_composite_afile.py [--data DIR] [--out data/a175531_composite.txt]
"""
import argparse, math, os, sys
from sympy import isprime

def typep(n, primes):
    return sorted(p for p in primes if n % (p * p - 1) == p)


DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')


def load(data):
    rows = {}  # n -> (primes, source)

    def add(n, primes, src):
        primes = sorted(primes)
        if n in rows:
            if rows[n][0] != primes:
                sys.exit(f'conflicting factorisations for {str(n)[:40]}...')
            return
        rows[n] = (primes, src)

    for name in ('cn2_nonrigid_multi.tsv', 'cn2_ktree_multi.tsv', 'cn2_divisible_5_7.tsv'):
        with open(os.path.join(data, name)) as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
            i_n, i_pr = hdr.index('n'), hdr.index('primes')
            for line in fh:
                f = line.rstrip('\n').split('\t')
                add(int(f[i_n]), [int(x) for x in f[i_pr].split()], name)
    return rows


def classify(n, primes):
    tp = typep(n, primes)
    if 5 in primes:
        return 1, tp
    if 7 in primes:
        return 2, tp
    return {4: 3, 3: 4, 2: 5}[len(tp)], tp


HEADER = """\
# Terms of A175531 (Carmichael numbers of order 2) of three special kinds: those divisible by 5 or 7,
# and those with k == p (mod p^2-1) at two, three or four primes p | k (and k == 1 (mod q^2-1) at
# every other prime q | k). Call such a p "type-p".
# These are EXAMPLES, not complete lists: each part comes from searches of particular moduli
# (complete searches of particular pools of primes in Parts 1 and 2; Howe's construction in Parts 3-5),
# not from a search of all k below a bound. The smallest term in each part is the smallest FOUND,
# not necessarily the smallest that exists. No term below 10^30 is of any of these kinds.
# Every number below was verified directly: its listed prime factors are distinct primes whose product
# is k, and for every prime p | k, k == 1 or k == p (mod p^2-1).
# Columns: k  log10(k)  type-p primes joined by '.' ('-' if none)  prime factors
#
# Part 1: {n1} terms divisible by 5. All are rigid (k == 1 (mod p^2-1) for every p | k).
#   The smallest, {d1} digits, was verified in Lean (CN2Lean/Counterexample5.lean).
# Part 2: {n2} terms divisible by 7. All are rigid.
#   The smallest, {d2} digits, was verified in Lean (CN2Lean/Counterexample7.lean).
# Part 3: {n3} term with four type-p primes, {d3} digits, {k3} prime factors.
# Part 4: {n4} terms with three type-p primes.
# Part 5: {n5} terms with two type-p primes.
# The terms of Parts 3-5 include all those with two or more type-p primes in the a-file of
# non-rigid examples already linked from A175531.
# Source and verifier: https://github.com/pudepiedj/carmichael-order2-nonrigid-multi
"""

TITLES = {1: 'Part 1: terms divisible by 5.',
          2: 'Part 2: terms divisible by 7.',
          3: 'Part 3: terms with four type-p primes.',
          4: 'Part 4: terms with three type-p primes.',
          5: 'Part 5: terms with two type-p primes.'}


def fmt(n, primes, tp):
    lg = f'{math.log10(n):.2f}' if n < 10 ** 300 else f'{len(str(n)) - 1 + math.log10(int(str(n)[:15]) / 10 ** 14):.2f}'
    return ' '.join([str(n), lg, '.'.join(map(str, tp)) or '-'] + [str(p) for p in primes])


def verify(path):
    part, counts, seen = None, {}, set()
    for line in open(path):
        if line.startswith('# Part ') and line.rstrip().endswith('.') and ':' in line and 'terms' in line and len(line) < 60:
            part = int(line.split()[2].rstrip(':'))
            continue
        if line.startswith('#') or not line.strip():
            continue
        f = line.split()
        n, tp, ps = int(f[0]), ([] if f[2] == '-' else [int(x) for x in f[2].split('.')]), [int(x) for x in f[3:]]
        assert n not in seen, 'duplicate'
        seen.add(n)
        assert len(set(ps)) == len(ps) and ps == sorted(ps), 'primes not distinct/sorted'
        assert all(isprime(p) for p in ps), 'non-prime factor'
        assert math.prod(ps) == n, 'product != k'
        assert n % 2 == 1 and len(ps) >= 2
        got = []
        for p in ps:
            r = n % (p * p - 1)
            assert r in (1, p), f'p={p} fails order-2 condition'
            if r == p:
                got.append(p)
        assert got == tp, 'type-p column wrong'
        assert abs(float(f[1]) - (len(str(n)) - 1 + math.log10(int(str(n)[:15]) / 10 ** 14))) < 0.006
        want = {1: lambda: 5 in ps, 2: lambda: 7 in ps and 5 not in ps, 3: lambda: len(tp) == 4,
                4: lambda: len(tp) == 3, 5: lambda: len(tp) == 2}[part]
        assert want(), f'row in wrong part {part}'
        assert n > 10 ** 30
        counts[part] = counts.get(part, 0) + 1
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default=DATA)
    ap.add_argument('--out', default=os.path.join(DATA, 'a175531_composite.txt'))
    a = ap.parse_args()
    rows = load(a.data)
    parts = {i: [] for i in range(1, 6)}
    for n, (primes, src) in rows.items():
        part, tp = classify(n, primes)
        parts[part].append((n, primes, tp))
    for v in parts.values():
        v.sort()
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    with open(a.out, 'w') as fh:
        fh.write(HEADER.format(n1=len(parts[1]), d1=len(str(parts[1][0][0])), n2=len(parts[2]),
                               d2=len(str(parts[2][0][0])), n3=len(parts[3]), d3=len(str(parts[3][0][0])),
                               k3=len(parts[3][0][1]), n4=len(parts[4]), n5=len(parts[5])))
        for i in range(1, 6):
            fh.write('#\n# ' + TITLES[i] + '\n')
            for n, primes, tp in parts[i]:
                fh.write(fmt(n, primes, tp) + '\n')
    src = {}
    for _, (_, s) in rows.items():
        src[s] = src.get(s, 0) + 1
    print('sources:', src)
    print('parts  :', {i: len(v) for i, v in parts.items()})
    print('verify :', verify(a.out))
    print('size   :', os.path.getsize(a.out), 'bytes ->', a.out)


if __name__ == '__main__':
    main()
