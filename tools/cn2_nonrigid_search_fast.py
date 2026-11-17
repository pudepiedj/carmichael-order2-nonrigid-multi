#!/usr/bin/env python3
"""cn2_nonrigid_search_fast.py -- search_multi() of cn2_nonrigid_search.py on the C engine cn2x/fast/cn2ss4.

    python3 cn2_nonrigid_search_fast.py L0 p0 [p1 ...] [--threads 32] [--rundir DIR] [--show 5] [--json runs.jsonl]
    python3 cn2_nonrigid_search_fast.py 33357220718421600 673 2017 2593        # the two triples of note 9.4
    python3 cn2_nonrigid_search_fast.py 2^7.3^3.5^2.7.11.13.17.19.29.31 1153   # Howe 1998: 53

Same problem, same verification and the same result dictionary as search_multi(): n = prod(T) * n0 with n0
squarefree over P(2, L0), n == 1 (mod q^2-1) for q | n0 and n == p (mod p^2-1) for p in T.  Only the
subset-product solver differs (see the header of cn2ss4.c): unequal quarters, hash join instead of sorting,
Barrett reduction.  The run is RESUMABLE: interrupt it (Ctrl-C, or touch RUNDIR/STOP) and repeat the same
command (any --threads); finished class-chunks are skipped.

There is no pool-size cap here beyond time and memory: cost ~ 2^(|P|/2) list entries, i.e. x1.41 per extra
pool prime.  --qmax (default 30) bounds the two big quarter tables at 2^qmax subsets (8 bytes per subset when built).
"""
import argparse
import itertools
import json
import math
import os
import subprocess
import sys
import time

from sympy import factorint, isprime
from sympy.ntheory.modular import crt

from cn2_systematic_fast_SS import calculate_P, parse_generator, phi_from_factors

HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.path.join(HERE, 'cn2x', 'fast', 'cn2ss4')


def howe_types(n, ps):
    out = {}
    for p in ps:
        r = n % (p * p - 1)
        out[p] = '1' if r == 1 else ('p' if r == p else None)
    return out


def setup(L0, L0f, p0s):
    """(L, Lf, pool, target, m) exactly as search_multi(); raises ValueError if T is not admissible on L0."""
    p0s = sorted(p0s)
    mods, res = [L0], [1]
    for p in p0s:
        if not isprime(p) or L0 % p == 0:
            raise ValueError(f"{p} is not prime or divides L0")
        mods.append(p * p - 1)
        res.append(p)
    for i in range(len(mods)):
        for j in range(i + 1, len(mods)):
            if (res[i] - res[j]) % math.gcd(mods[i], mods[j]):
                raise ValueError(f"incompatible: moduli {mods[i]}, {mods[j]}")
    for p in p0s:
        if any((q * q - 1) % p == 0 for q in p0s if q != p):
            raise ValueError(f"{p} divides another type-p prime's p^2-1")
    L = math.lcm(*mods)
    Lf = factorint(L)
    P = [q for q in calculate_P(L0f) if all((p * p - 1) % q for p in p0s) and q not in p0s]
    c = int(crt(mods, res)[0])
    m = math.prod(p0s)
    return L, Lf, P, c * pow(m, -1, L) % L, m


def plan(Lf, n, qmax=30, list_bits=16, max_A_bits=27, yield_bits=4):
    """Quarter sizes and the coprime split L = A * G * E_1 * ...

    A: class modulus, phi(A) as close as possible to 2^(n/2 - list_bits) so the per-class build list is
       ~2^list_bits entries (its hash table then stays in L2).  The engine uses the true |H| = |<P mod A>|.
    quarters: the big ones (Q2, Q4) as large as qmax allows, so a bucket lookup yields a long run."""
    half1 = n // 2
    half2 = n - half1
    want = max(0, half1 - list_bits)
    comps = sorted(p ** e for p, e in Lf.items())
    phis = {p ** e: (p - 1) * p ** (e - 1) for p, e in Lf.items()}
    best = None
    small = [c for c in comps if c < (1 << max_A_bits)]
    for r in range(0, min(len(small), 6) + 1):
        for sub in itertools.combinations(small, r):
            A = math.prod(sub)
            if A >= (1 << max_A_bits):
                continue
            hb = math.log2(math.prod(phis[c] for c in sub)) if sub else 0.0
            key = (abs(hb - want), A)
            if best is None or key < best[0]:
                best = (key, A, sub, hb)
    _, A, sub, hb = best
    rest = sorted((c for c in comps if c not in sub), reverse=True)
    gbest = (0, 1, ())                            # join modulus G < 2^32 with the largest phi: fewest false candidates
    for r in range(1, min(len(rest), 7) + 1):
        for gs in itertools.combinations(rest, r):
            Gm = math.prod(gs)
            if Gm < (1 << 32):
                gbest = max(gbest, (math.prod(phis[c] for c in gs), Gm, gs))
    groups = [gbest[1]]
    for c in rest:                                # the remainder: first-fit into moduli < 2^62, checked per candidate
        if c in gbest[2]:
            continue
        if c >= (1 << 62):
            raise ValueError(f"prime power {c} too large for the engine")
        for g in range(1, len(groups)):
            if groups[g] * c < (1 << 62):
                groups[g] *= c
                break
        else:
            groups.append(c)

    def big(half):                                # size of the bucketed quarter
        q = max(half // 2, math.ceil(hb) + yield_bits)
        return max(1, min(qmax, half - 1, q))
    q2, q4 = big(half1), big(half2)
    return dict(A=A, G=groups[0], E=groups[1:], qs=[half1 - q2, q2, half2 - q4, q4], h_bits=hb)


def build_engine():
    src = ENGINE + '.c'
    if not os.path.exists(ENGINE) or os.path.getmtime(ENGINE) < os.path.getmtime(src):
        subprocess.run(['clang', '-O3', '-mcpu=native', '-o', ENGINE, src, '-lpthread', '-lm'], check=True)


def search_multi_fast(L0, L0f, p0s, rundir=None, threads=32, verbose=True, qmax=30, list_bits=16, progress=30,
                      yield_bits=4):
    p0s = sorted(p0s)
    L, Lf, P, target, m = setup(L0, L0f, p0s)
    n = len(P)
    predicted = 2.0 ** n / phi_from_factors(Lf)
    if n < 4:
        raise ValueError(f"pool of {n} primes is too small")
    pl = plan(Lf, n, qmax, list_bits, yield_bits=yield_bits)
    A, G, E, qs = pl['A'], pl['G'], pl['E'], pl['qs']
    assert A * G * math.prod(E) == L
    job = f"{n} {' '.join(map(str, qs))}\n{' '.join(map(str, P))}\n{A} {target % A}\n{G} {target % G}\n{len(E)}\n"
    job += ''.join(f"{e} {target % e}\n" for e in E)
    if rundir is None:
        rundir = os.path.join(HERE, 'cn2x', 'runs', f"ss4_{L0}_{'_'.join(map(str, p0s))}")
    os.makedirs(rundir, exist_ok=True)
    jp = os.path.join(rundir, 'job.txt')
    if os.path.exists(jp):
        if open(jp).read() != job:
            raise RuntimeError(f"{jp} holds a different job; use another --rundir or remove that directory")
    else:
        with open(jp, 'w') as fh:
            fh.write(job)
    build_engine()
    if verbose:
        print(f"L = {L} ({L.bit_length()} bits)  |P| = {n}  predicted {predicted:.3f}  "
              f"work ~2^{math.log2(2.0 ** (qs[0] + qs[1]) + 2.0 ** (qs[2] + qs[3])):.1f} list entries\nrundir {rundir}",
              flush=True)
    t0 = time.time()
    rc = subprocess.run([ENGINE, rundir, str(threads), str(progress)],
                        stderr=None if verbose else subprocess.DEVNULL).returncode
    seconds = time.time() - t0
    if rc not in (0, 2):
        raise RuntimeError(f"cn2ss4 failed with exit code {rc}")
    offs = [0, qs[0], qs[0] + qs[1], qs[0] + qs[1] + qs[2]]
    values = set()
    op = os.path.join(rundir, 'out.txt')
    for line in open(op) if os.path.exists(op) else []:
        f = line.split()
        if len(f) == 5 and f[0] == 'H':
            n0 = 1
            for q in range(4):
                mask = int(f[q + 1])
                for b in range(qs[q]):
                    if mask >> b & 1:
                        n0 *= P[offs[q] + b]
            values.add(n0)
    found = []
    for n0 in sorted(values):
        N = m * n0
        ps = sorted([q for q in P if n0 % q == 0] + p0s)
        ty = howe_types(N, ps)
        ok = math.prod(ps) == N and all(ty.values()) and [p for p in ps if ty[p] == 'p'] == p0s
        found.append((N, ps, ok))
    found.sort()
    return dict(L0=L0, p0s=p0s, L=L, P=n, predicted=predicted, count=len(found), collected=len(found),
                verified=sum(ok for _, _, ok in found), seconds=seconds, found=found, complete=(rc == 0),
                rundir=rundir, plan=pl)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('L0', help="generator: 33357220718421600 or 2^5.3^4.5^2...")
    ap.add_argument('p0', type=int, nargs='+', help="the type-p primes T (one or more)")
    ap.add_argument('--threads', type=int, default=32)
    ap.add_argument('--rundir')
    ap.add_argument('--qmax', type=int, default=30, help='log2 of the largest bucketed quarter table')
    ap.add_argument('--list-bits', type=int, default=16, help='log2 of the per-class build list aimed for')
    ap.add_argument('--yield-bits', type=int, default=4, help='log2 of the entries aimed for per bucket lookup')
    ap.add_argument('--show', type=int, default=5)
    ap.add_argument('--json', help='append one JSON line')
    a = ap.parse_args()
    L0, L0f = parse_generator(a.L0)
    r = search_multi_fast(L0, L0f, a.p0, a.rundir, a.threads, True, a.qmax, a.list_bits, yield_bits=a.yield_bits)
    state = '' if r['complete'] else '  INCOMPLETE -- rerun the same command to resume'
    print(f"L0 = {L0}  T = {r['p0s']}  |P| = {r['P']}  predicted {r['predicted']:.2f}  found {r['count']}  "
          f"verified {r['verified']}  {r['seconds']:.1f}s{state}")
    for N, ps, ok in r['found'][:a.show]:
        print(f"  {'OK ' if ok else 'BAD'} 10^{math.log10(N):.2f}  k={len(ps)}  n = {N}")
        print(f"      {', '.join(map(str, ps))}")
    if a.json and r['complete']:
        with open(a.json, 'a') as fh:
            rec = {k: v for k, v in r.items() if k not in ('found', 'plan')}
            rec['values'] = [str(N) for N, _, _ in r['found']]
            fh.write(json.dumps(rec) + '\n')


if __name__ == '__main__':
    main()
