#!/usr/bin/env python3
"""cn2_rich_multi_scout.py -- can a very fertile generator carry k type-p primes and still be searched in part?

    python3 cn2_rich_multi_scout.py [--p0max 200000] [--top 30] [--kmax 6]

For a generator L0 and a set T of type-p primes (Howe's construction, as in cn2_multi_scout.py):
    P' = pool of L0 minus primes dividing some p^2-1 or in T,   beta = log2 phi(lcm(L0, p^2-1 : p in T)),
    F_T = |P'| - beta   (log2 of the expected count E).
Cost to an expected first hit, in list entries:
    exhaustive engine (cn2ss4), or any partial run of it:  2^(|P'|/2 - F_T)
    Wagner k-tree with 2^m lists (lean_support/ktree.c is m = 2), one sub-pool of 2^m * beta/(m+1) primes:
        about 2^(beta/(m+1)) entries per list for about one hit per run, if |P'| >= 2^m * beta/(m+1).
So a partial run only beats the full one when F_T is large, and a very large pool makes beta, not |P'|, the cost.
The generators are hand-picked rich ones (L_5 of CN2_Divisible_by_5_or_7.md and larger), not the box.
"""
import argparse, itertools, math, os, sys, time
from sympy import factorint, primerange

from cn2_nonrigid_scout import log2phi, merge

HERE = os.path.dirname(os.path.abspath(__file__))


def pool(f):
    """Howe's pool of L = prod p^e, by cn2x/fast/pool2 (threaded C, no divisor list; d = q-1 < 2^62);
    same result as cn2_nonrigid_scout.pool, which stores every divisor and cannot reach these generators"""
    import subprocess
    d = os.path.join(HERE, 'cn2x', 'fast')
    exe = os.path.join(d, 'pool2')
    if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(exe + '.c'):
        subprocess.run(['clang', '-O3', '-mcpu=native', '-o', exe, exe + '.c', '-lpthread'], check=True)
    args = [str(x) for p, e in sorted(f.items()) for x in (p, e)]
    out = subprocess.run([exe, str(os.cpu_count() or 8), *args], capture_output=True, text=True, check=True).stdout
    return sorted(int(x) for x in out.split())


def v2(x):
    return (x & -x).bit_length() - 1


def gens():
    def prod_primes(lo, hi):
        return {p: 1 for p in primerange(lo, hi + 1)}
    L5 = {2: 9, 3: 5, 7: 2, **prod_primes(11, 61)}
    yield 'L5 (5-free)', L5
    yield 'L5*5^2', {**L5, 5: 2}
    yield '2^9 3^5 5^2 7^2 11..67', {2: 9, 3: 5, 5: 2, 7: 2, **prod_primes(11, 67)}
    yield '2^10 3^5 5^2 7^2 11..71', {2: 10, 3: 5, 5: 2, 7: 2, **prod_primes(11, 71)}
    yield '2^10 3^6 5^2 7^2 11^2 13..73', {2: 10, 3: 6, 5: 2, 7: 2, 11: 2, **prod_primes(13, 73)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--p0max', type=int, default=200000)
    ap.add_argument('--top', type=int, default=30)
    ap.add_argument('--kmax', type=int, default=5)
    ap.add_argument('--gen', action='append', help='generator as 2^10.3^6.5^2... (repeatable; default: the built-in list)')
    a = ap.parse_args()
    if a.gen:
        global gens
        specs = [(g, factorint(int(g)) if g.isdigit() else
                  {int(t.partition('^')[0]): int(t.partition('^')[2] or 1) for t in g.split('.')}) for g in a.gen]
        gens = lambda: iter(specs)
    cand = [p for p in primerange(3, a.p0max) if p % 24 == 1]
    fac = {}
    for p in cand:
        fm, fp = factorint(p - 1), factorint(p + 1)
        fac[p] = merge(fm, fp)
        fac[p][2] = fm.get(2, 0) + fp.get(2, 0)
    for name, f in gens():
        t0 = time.time()
        L0 = math.prod(p ** e for p, e in f.items())
        P = pool(f)
        lp = log2phi(f)
        print(f"\n{name}: L0 ~ 2^{math.log2(L0):.1f}, pool {len(P)}, beta {lp:.1f}, F {len(P) - lp:+.1f}, "
              f"P/4 = {len(P) / 4:.1f}  ({time.time() - t0:.0f}s)", flush=True)
        adm = [p for p in cand if L0 % p and (p - 1) % math.gcd(L0, p * p - 1) == 0]
        cost1 = {p: log2phi(merge(f, fac[p])) - lp + sum(1 for q in P if (p * p - 1) % q == 0 or q == p) for p in adm}
        byv = {}
        for p in sorted(adm, key=cost1.get):
            byv.setdefault(v2(p - 1), []).append(p)
        best = {}
        kill = {p: sum(1 << i for i, q in enumerate(P) if (p * p - 1) % q == 0 or q == p) for p in adm}
        for v, ps in byv.items():
            ps = sorted(ps[:a.top])
            ok = {(x, y): (y - x) % math.gcd(x * x - 1, y * y - 1) == 0 and (x * x - 1) % y and (y * y - 1) % x
                  for x, y in itertools.combinations(ps, 2)}
            for k in range(1, a.kmax + 1):
                for T in itertools.combinations(ps, k):
                    if not all(ok[pr] for pr in itertools.combinations(T, 2)):
                        continue
                    fl = f
                    for p in T:
                        fl = merge(fl, fac[p])
                    m = 0
                    for p in T:
                        m |= kill[p]
                    Pe = len(P) - bin(m).count('1')
                    beta = log2phi(fl)
                    if k not in best or beta < best[k][1]:
                        best[k] = (Pe - beta, beta, Pe, T)
        for k in sorted(best):
            FT, beta, Pe, T = best[k]
            ex = max(Pe / 2 - FT, 0)
            wg = [(m, beta / (m + 1), 2 ** m * beta / (m + 1)) for m in range(2, 8) if 2 ** m * beta / (m + 1) <= Pe]
            w = f"2^m={2 ** wg[0][0]} lists of 2^{wg[0][1]:.1f}" if wg else "no k-tree fits"
            w2 = ", ".join(f"{2 ** m}:2^{c:.1f}" for m, c, _ in wg)
            print(f"  k={k}: beta {beta:6.1f}  F_T {FT:+7.1f}  pool {Pe}  exhaustive/partial 2^{ex:.1f};  "
                  f"k-tree {w}  [{w2}]  T={list(T)}", flush=True)


if __name__ == '__main__':
    main()
