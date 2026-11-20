#!/usr/bin/env python3
"""cn2_ktree.py -- Wagner's k-tree (generalised birthday) search for Howe-construction order-2 Carmichael numbers.

    python3 cn2_ktree.py L0 p0 [p1 ...] [--m 4] [--b 26] [--runs 4] [--seed 1] [--plan] [--json out.jsonl]
    python3 cn2_ktree.py 2^10.3^6.5^2.7^2.11^2.13.17.19.23.29.31.37.41.43.47.53.59.61.67.71.73 15361 --m 4 --b 27

Problem (as search_multi): n = prod(T) * n0, n0 squarefree over the pool P' of L0, with
prod(S) == target (mod L), L = lcm(L0, p^2-1 : p in T).  Every hit is rebuilt and checked exactly
(n == 1 mod q^2-1 at each pool prime, n == p mod p^2-1 at each p in T, type-p set == T).

Method (Wagner, CRYPTO 2002, multiplicative form; lean_support/ktree.c is the case m = 2):
  * L is split into coprime parts M_1 .. M_m (join moduli) and filter moduli F_1 .. (the rest), all < 2^32.
  * 2^m leaf lists, each all 2^b subsets of b distinct pool primes (random choice per run), leaf i scaled by a
    random unit u_i with prod u_i == target^-1 (mod L).
  * Level l merges pairs of lists keeping a*b == 1 (mod M_l); the last merge is filtered by the F's.
Each run finds only a slice of the solutions (it can never certify absence).  Expected hits per run:
    E_run = 2^(2^m b) / phi(L) / prod_l phi(M_l)^(2^(m-l) - 1).
List sizes stay near 2^b when phi(M_l) ~ 2^b; memory ~ (2^(m+1)) * 2^b * 4 bytes of provenance plus residues.
This is a numpy prototype: practical to b ~ 28.
"""
import argparse, json, math, os, random, sys, time

import numpy as np
from sympy import factorint, isprime
from sympy.ntheory.modular import crt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cn2_rich_multi_scout import pool                    # noqa: E402  (C pool, cn2x/fast/pool2)
from cn2_nonrigid_search_fast import howe_types         # noqa: E402

LIM = 1 << 32


def parse_L(s):
    if s.isdigit():
        return int(s)
    n = 1
    for t in s.split('.'):
        b, _, e = t.partition('^')
        n *= int(b) ** int(e or 1)
    return n


def phi_bits(pp):
    (p, e), = factorint(pp).items()
    return math.log2(p - 1) + (e - 1) * math.log2(p)


def setup(L0, p0s):
    """L, prime-power list of L, pool P', target (prod S == target mod L), prod(T); as cn2_nonrigid_search_fast.setup"""
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
    P = [q for q in pool(factorint(L0)) if all((p * p - 1) % q for p in p0s) and q not in p0s]
    c = int(crt(mods, res)[0])
    mT = math.prod(p0s)
    return L, Lf, P, c * pow(mT, -1, L) % L, mT


def plan(Lf, m, b, bR=None, lim=None):
    """join moduli M_1..M_m with phi ~ the list size each, the rest as filter moduli; all < lim (2^32 numpy, 2^40 C).
    Left leaves have b primes, right leaves bR (default b); lists are kept near 2^((b + bR)/2)."""
    LIM = lim or globals()['LIM']
    bR = b if bR is None else bR
    items = sorted(((p ** e, phi_bits(p ** e)) for p, e in Lf.items()), key=lambda t: -t[1])
    if items[0][0] >= LIM:
        raise ValueError(f"prime power {items[0][0]} of L is >= {LIM}")
    joins, s = [], (b + bR) / 2
    rng = random.Random(12345)                    # deterministic: the same plan on every rerun
    for l in range(m):
        tgt = s if l < m - 1 else math.log2(LIM) - 0.5   # keep list sizes level; the last join takes all it can
        best = None
        for trial in range(4000):                 # randomised first-fit, closest to the target from below wins
            order = items if trial == 0 else rng.sample(items, len(items))
            M, bits, used = 1, 0.0, []
            for it in order:
                if bits + it[1] <= tgt + 0.1 and M * it[0] < LIM:
                    M *= it[0]; bits += it[1]; used.append(it)
            if best is None or bits > best[1]:
                best = (M, bits, used)
            if tgt - bits < 0.05:
                break
        M, bits, used = best
        for it in used:
            items.remove(it)
        joins.append((M, bits))
        s = 2 * s - bits
    filters = []
    for pp, bits in items:
        for i, (M, fb) in enumerate(filters):
            if M * pp < LIM:
                filters[i] = (M * pp, fb + bits)
                break
        else:
            filters.append((pp, bits))
    return joins, filters


def expected(m, b, joins, beta, bR=None):
    bR = b if bR is None else bR
    lg = 2 ** (m - 1) * (b + bR) - beta - sum((2 ** (m - l) - 1) * jb for l, (_, jb) in enumerate(joins, 1))
    sizes, s = [], (b + bR) / 2
    for _, jb in joins:
        s = 2 * s - jb
        sizes.append(s)
    return lg, sizes


def leaf(primes, mods, scale):
    """all subset products of primes: forward and inverse residues per modulus, and the subset mask"""
    k = len(primes)
    fw, iv = [], []
    for M, u in zip(mods, scale):
        a = np.array([u % M], dtype=np.uint64)
        ai = np.array([pow(u, -1, M)], dtype=np.uint64)
        for p in primes:
            a = np.concatenate([a, a * np.uint64(p % M) % np.uint64(M)])
            ai = np.concatenate([ai, ai * np.uint64(pow(p, -1, M)) % np.uint64(M)])
        fw.append(a.astype(np.uint32)); iv.append(ai.astype(np.uint32))
    mask = np.zeros(1, dtype=np.uint32)
    for i in range(k):
        mask = np.concatenate([mask, mask | np.uint32(1 << i)])
    return fw, iv, mask


_LIB = None
THREADS = os.cpu_count() or 8


def lib():
    """build (if stale) and load cn2x/fast/libktjoin.dylib"""
    global _LIB
    if _LIB is None:
        import ctypes, subprocess
        d = os.path.join(HERE, 'cn2x', 'fast')
        src, so = os.path.join(d, 'ktjoin.c'), os.path.join(d, 'libktjoin.dylib')
        if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
            subprocess.run(['clang', '-O3', '-mcpu=native', '-shared', '-fPIC', '-o', so, src, '-lpthread'], check=True)
        _LIB = ctypes.CDLL(so)
        u32 = np.ctypeslib.ndpointer(np.uint32, flags='C_CONTIGUOUS')
        _LIB.kt_join.restype = ctypes.c_int64
        _LIB.kt_join.argtypes = [u32, ctypes.c_int64, u32, ctypes.c_int64, u32, u32, ctypes.c_int64, ctypes.c_int]
        _LIB.kt_mulgather.restype = None
        _LIB.kt_mulgather.argtypes = [u32, u32, u32, u32, ctypes.c_int64, ctypes.c_uint32, u32, ctypes.c_int]
    return _LIB


def join(A, B, mods, j):
    """pairs (a, b) with a*b == 1 mod mods[j]; returns merged node keeping moduli j+1.. and provenance"""
    fa, ia = A['fw'], A['iv']
    fb, ib = B['fw'], B['iv']
    ka, kb = ia[0], fb[0]                     # a^-1 == b  (mod M_j); index 0 is always the current join modulus
    L = lib()
    cap = max(1 << 16, int(2.5 * len(ka) * len(kb) / mods[j]) + (1 << 16))
    while True:
        xa = np.empty(cap, np.uint32); xb = np.empty(cap, np.uint32)
        tot = L.kt_join(ka, len(ka), kb, len(kb), xa, xb, cap, THREADS)
        if tot <= cap:
            break
        cap = tot
    xa, xb = xa[:tot].copy(), xb[:tot].copy()
    fw, iv = [], []
    for t in range(1, len(fa)):
        M = mods[j + t]
        o1 = np.empty(tot, np.uint32); o2 = np.empty(tot, np.uint32)
        L.kt_mulgather(fa[t], fb[t], xa, xb, tot, M, o1, THREADS)
        L.kt_mulgather(ia[t], ib[t], xa, xb, tot, M, o2, THREADS)
        fw.append(o1); iv.append(o2)
    return dict(fw=fw, iv=iv, pa=xa, pb=xb, A=A, B=B)


def strip(node):
    """drop residues of a finished child (provenance kept)"""
    node['fw'] = node['iv'] = None


def subsets(node, idx):
    if 'mask' in node:
        m = int(node['mask'][idx])
        return [q for i, q in enumerate(node['primes']) if m >> i & 1]
    return subsets(node['A'], int(node['pa'][idx])) + subsets(node['B'], int(node['pb'][idx]))


def run(L, Lf, P, target, m, b, rng, verbose=True):
    joins, filters = plan(Lf, m, b)
    mods = [M for M, _ in joins] + [M for M, _ in filters]
    if math.prod(mods) != L:
        raise RuntimeError("partition of L is wrong")
    need = (2 ** m) * b
    if need > len(P):
        raise ValueError(f"2^m * b = {need} primes needed, pool has {len(P)}")
    prs = rng.sample(P, need)
    us = []
    for _ in range(2 ** m - 1):
        while True:
            u = rng.randrange(2, L)
            if math.gcd(u, L) == 1:
                break
        us.append(u)
    us.append(pow(target * math.prod(us), -1, L))       # prod(u_i) == target^-1, so prod(u_i) * prod(S) == 1 <=> prod(S) == target
    t0 = time.time()

    def build(level, i):
        """depth-first: only one branch's residues are alive at a time (provenance of every node is kept)"""
        if level == 0:
            ps = sorted(prs[i * b:(i + 1) * b])
            fw, iv, mask = leaf(ps, mods, [us[i]] * len(mods))
            return dict(fw=fw, iv=iv, mask=mask, primes=ps)
        A = build(level - 1, 2 * i)
        B = build(level - 1, 2 * i + 1)
        node = join(A, B, mods, level - 1)
        strip(A); strip(B)
        if verbose and i == 0:
            print(f"    level {level}: first list 2^{math.log2(max(1, len(node['pa']))):.1f}  ({time.time() - t0:.0f}s)",
                  flush=True)
        return node

    root = build(m, 0)
    ok = np.ones(len(root['pa']), dtype=bool)
    for t in range(len(filters)):
        ok &= root['fw'][t] == 1
    return [subsets(root, int(i)) for i in np.nonzero(ok)[0]]


def engine_c():
    import subprocess
    d = os.path.join(HERE, 'cn2x', 'fast')
    exe = os.path.join(d, 'cn2kt')
    if not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(exe + '.c'):
        subprocess.run(['clang', '-O3', '-mcpu=native', '-o', exe, exe + '.c', '-lpthread'], check=True)
    return exe


def footprint(m, b, bR, joins, filters):
    """approximate peak memory and prov disk (bytes) of cn2kt for this plan"""
    W = [(M - 1).bit_length() for M, _ in joins] + [(M - 1).bit_length() for M, _ in filters]
    R = [math.ceil(sum(W[l:]) / 8) for l in range(len(W) + 1)]
    _, sizes = expected(m, b, joins, 0.0, bR)
    N = [2 ** x for x in sizes]                                  # N[l-1] = entries at level l
    peak = sum(N[l - 1] * R[l] for l in range(1, m)) + (N[0] * R[1] if m > 1 else 0) \
        + (N[1] * (R[2] + 10) if m > 2 else 0) + 2 ** bR * math.ceil((W[0] + 32) / 8)
    disk = sum(2 ** (m - l) * N[l - 1] * 10 for l in range(2, m))
    return peak, disk


def run_c(L, P, target, m, b, bR, joins, filters, rng, rundir, threads=32):
    """one run on cn2x/fast/cn2kt; returns the hit subsets (lists of pool primes)"""
    import subprocess
    mods = [M for M, _ in joins] + [M for M, _ in filters]
    if math.prod(mods) != L:
        raise RuntimeError("partition of L is wrong")
    nl = 2 ** m
    bs = [b if i % 2 == 0 else bR for i in range(nl)]
    if sum(bs) > len(P):
        raise ValueError(f"{sum(bs)} primes needed, pool has {len(P)}")
    prs = rng.sample(P, sum(bs))
    groups, o = [], 0
    for x in bs:
        groups.append(sorted(prs[o:o + x])); o += x
    us = []
    for _ in range(nl - 1):
        while True:
            u = rng.randrange(2, L)
            if math.gcd(u, L) == 1:
                break
        us.append(u)
    us.append(pow(target * math.prod(us), -1, L))
    os.makedirs(rundir, exist_ok=True)
    with open(os.path.join(rundir, 'job.txt'), 'w') as fh:
        fh.write(f"{len(mods)} {m} 2\n{' '.join(map(str, mods))}\n")
        for g in groups:
            fh.write(f"{len(g)} {' '.join(map(str, g))}\n")
        for u in us:
            fh.write(' '.join(str(u % M) for M in mods) + '\n')
    hp = os.path.join(rundir, 'hits.txt')
    if os.path.exists(hp):
        os.remove(hp)
    rc = subprocess.run([engine_c(), rundir, str(threads)]).returncode
    if rc != 0:
        raise RuntimeError(f"cn2kt failed with exit code {rc}")
    hits = []
    for line in open(hp):
        f = line.split()
        if f and f[0] == 'H':
            S = []
            for g, mk in zip(groups, map(int, f[1:])):
                S += [q for t, q in enumerate(g) if mk >> t & 1]
            hits.append(S)
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('L0')
    ap.add_argument('p0s', type=int, nargs='+')
    ap.add_argument('--m', type=int, default=4, help='2^m leaf lists')
    ap.add_argument('--b', type=int, default=24, help='primes per leaf list (2^b entries)')
    ap.add_argument('--runs', type=int, default=1)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--plan', action='store_true', help='show the split and expected yield only')
    ap.add_argument('--json', help='append one record per run')
    ap.add_argument('--engine', choices=['numpy', 'c'], default='numpy',
                    help='c = cn2x/fast/cn2kt: bit-packed, provenance on disk, moduli < 2^40 (for lists of 2^30+)')
    ap.add_argument('--bR', type=int, help='primes per RIGHT leaf (stored); default --b. The left leaves are never stored')
    ap.add_argument('--threads', type=int, default=32)
    ap.add_argument('--rundir', help='work directory for --engine c (provenance files; default under cn2x/runs)')
    a = ap.parse_args()
    if a.engine == 'numpy' and a.bR not in (None, a.b):
        raise SystemExit('--bR needs --engine c')
    L0 = parse_L(a.L0)
    t0 = time.time()
    L, Lf, P, target, mT = setup(L0, a.p0s)
    beta = sum(math.log2(p - 1) + (e - 1) * math.log2(p) for p, e in Lf.items())
    bR = a.b if a.bR is None else a.bR
    joins, filters = plan(Lf, a.m, a.b, bR, lim=(1 << 40) if a.engine == 'c' else LIM)
    lg, sizes = expected(a.m, a.b, joins, beta, bR)
    print(f"L0 = {L0}  T = {sorted(a.p0s)}  pool {len(P)} (setup {time.time() - t0:.0f}s)  beta {beta:.1f}  "
          f"F_T {len(P) - beta:+.1f}\n  {2 ** a.m} lists of {a.b}/{bR} primes ({2 ** (a.m - 1) * (a.b + bR)} used); join phi bits "
          f"{[round(jb, 1) for _, jb in joins]}; filter phi bits {[round(fb, 1) for _, fb in filters]}\n"
          f"  list sizes by level 2^{[round(s, 1) for s in sizes]}; expected hits per run 2^{lg:.2f} = {2 ** lg:.3g}",
          flush=True)
    if a.engine == 'c':
        pk, dk = footprint(a.m, a.b, bR, joins, filters)
        print(f"  cn2kt: peak memory ~{pk / 1e9:.0f} GB, provenance on disk ~{dk / 1e9:.0f} GB", flush=True)
    if a.plan:
        return
    rng = random.Random(a.seed)
    allv = set()
    for r in range(a.runs):
        t1 = time.time()
        if a.engine == 'c':
            rd = a.rundir or os.path.join(HERE, 'cn2x', 'runs', f"kt_{L0}_{'_'.join(map(str, sorted(a.p0s)))}")
            hits = run_c(L, P, target, a.m, a.b, bR, joins, filters, rng, rd, a.threads)
        else:
            hits = run(L, Lf, P, target, a.m, a.b, rng)
        vals = []
        for S in hits:
            N = mT * math.prod(S)
            ps = sorted(S + sorted(a.p0s))
            ty = howe_types(N, ps)
            good = len(set(ps)) == len(ps) and all(ty.values()) and [p for p in ps if ty[p] == 'p'] == sorted(a.p0s)
            vals.append((N, good))
        nv = sum(g for _, g in vals)
        allv |= {N for N, g in vals if g}
        print(f"  run {r + 1}: {len(vals)} hits ({nv} verified) in {time.time() - t1:.0f}s; expected {2 ** lg:.3g}"
              + (f"; smallest 10^{math.log10(min(N for N, _ in vals)):.1f} with "
                 f"{min((N, len(S)) for (N, _), S in zip(vals, hits))[1] + len(a.p0s)} primes" if vals else ""),
              flush=True)
        if nv != len(vals):
            raise RuntimeError("a hit failed exact verification")
        if a.json:
            with open(a.json, 'a') as fh:
                fh.write(json.dumps(dict(L0=str(L0), p0s=sorted(a.p0s), m=a.m, b=a.b, bR=bR, engine=a.engine, seed=a.seed, run=r + 1,
                                         expected=2 ** lg, count=len(vals), verified=nv,
                                         seconds=time.time() - t1, values=[str(N) for N, _ in vals])) + '\n')
    print(f"{a.runs} run(s): {len(allv)} distinct verified numbers")


if __name__ == '__main__':
    main()
