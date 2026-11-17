#!/usr/bin/env python3
"""cn2_multi_scout.py -- price compatible SETS of k type-p primes on the box generators (k = 3, 4, ...).

    python3 cn2_multi_scout.py --k 4 [--fmin 33] [--pmax 110] [--p0max 30000] [--top 20] [--out cn2_quad_sets.json]

Generalises cn2_twin_scout.py.  For a generator L0 and a set T of type-p primes the count model is
    E = 2^|P'| / phi(lcm(L0, p^2-1 : p in T)),      P' = P(2, L0) minus primes dividing some p^2-1 or in T,
and the search cost on cn2ss4 is ~2^(|P'|/2) list entries.  Per generator only the --top cheapest admissible
single primes are combined (a set is never cheaper than its members), so this is a scout, not a census.
Reports the sets ranked by E, by smallest pool with E >= 1, and by E per unit of work.
"""
import argparse, csv, itertools, json, math, os
from sympy import factorint, primerange

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--k', type=int, default=4)
    ap.add_argument('--fmin', type=float, default=33.0)
    ap.add_argument('--pmax', type=int, default=110)
    ap.add_argument('--p0max', type=int, default=30000)
    ap.add_argument('--top', type=int, default=20)
    ap.add_argument('--emin', type=float, default=0.01)
    ap.add_argument('--out')
    a = ap.parse_args()
    from cn2_nonrigid_scout import pool, log2phi, merge
    rows = [r for r in csv.DictReader(open(os.path.join(HERE, 'cn2_box_fecundity.csv')))
            if int(r['P']) <= a.pmax and float(r['F']) >= a.fmin]
    cand = [p for p in primerange(3, a.p0max) if p % 24 == 1]
    fac = {}
    for p in cand:
        fm, fp = factorint(p - 1), factorint(p + 1)
        fac[p] = merge(fm, fp)
        fac[p][2] = fm.get(2, 0) + fp.get(2, 0)
    v2 = {p: factorint(p - 1)[2] for p in cand}
    print(f"{len(rows):,} generators with F >= {a.fmin} and |P| <= {a.pmax}; {len(cand):,} candidate primes < {a.p0max}; "
          f"k = {a.k}", flush=True)
    out = []
    for n, r in enumerate(rows):
        L0 = int(r['L']); f = factorint(L0); P = pool(f); lp = log2phi(f); F = float(r['F'])
        adm = [p for p in cand if L0 % p and (p - 1) % math.gcd(L0, p * p - 1) == 0]
        d1 = {p: log2phi(merge(f, fac[p])) - lp + sum(1 for q in P if (p * p - 1) % q == 0 or q == p) for p in adm}
        byv = {}
        for p in sorted(adm, key=d1.get):
            byv.setdefault(v2[p], []).append(p)
        for v, ps in byv.items():
            ps = sorted(ps[:a.top])
            ok = {(x, y): (y - x) % math.gcd(x * x - 1, y * y - 1) == 0 and (x * x - 1) % y and (y * y - 1) % x
                  for x, y in itertools.combinations(ps, 2)}
            for T in itertools.combinations(ps, a.k):
                if not all(ok[pr] for pr in itertools.combinations(T, 2)):
                    continue
                fl = f
                for p in T:
                    fl = merge(fl, fac[p])
                Pe = sum(1 for q in P if q not in T and all((p * p - 1) % q for p in T))
                E = 2.0 ** (Pe - log2phi(fl))
                if E >= a.emin:
                    out.append(dict(L0=str(L0), T=list(T), P=Pe, E=E, F=F))
        if (n + 1) % 200 == 0:
            print(f"  [{n + 1:,}/{len(rows):,}] kept {len(out):,}", flush=True)
    print(f"{len(out):,} sets with E >= {a.emin}; total E {sum(d['E'] for d in out):.1f}")

    def show(title, ds):
        print(title)
        for d in ds:
            print(f"  E={d['E']:9.3f}  |P'|={d['P']:>3}  work 2^{d['P'] / 2:.1f}  F={d['F']:6.2f}  T={d['T']}  L0={d['L0']}")
    show("largest E:", sorted(out, key=lambda d: -d['E'])[:10])
    show("smallest pool with E >= 1:", sorted([d for d in out if d['E'] >= 1], key=lambda d: (d['P'], -d['E']))[:15])
    show("smallest pool with E >= 0.25:", sorted([d for d in out if d['E'] >= 0.25], key=lambda d: (d['P'], -d['E']))[:10])
    show("best E per unit work (E / 2^(|P'|/2)):", sorted(out, key=lambda d: -(math.log2(d['E']) - d['P'] / 2))[:10])
    if a.out:
        out.sort(key=lambda d: -d['E'])
        json.dump(out, open(a.out, 'w'), indent=0)
        print(f"-> {a.out}")


if __name__ == '__main__':
    main()
