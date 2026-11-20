#!/usr/bin/env python3
"""
CN2 systematic generator search: chunked, 64-bit-safe MITM (17 November 2026).

The same search as cn2_systematic_fast.py, but with fixes for three problems
that appear once the pool exceeds about 64 primes or L exceeds 2^32.

1. Memory. The original stores both halves (2^34 entries each for |P| = 68)
   plus a sort permutation, a sorted copy and two per-entry count arrays:
   about 700 GiB peak. Here only the stored half is kept, as one sorted
   uint64 array (2^STORED_BITS entries; 32 GiB at the default 32). The other
   half is streamed in chunks of 2^LOW_BITS, so it costs no memory.
2. Index overflow. The original keeps subset indices in uint32, which wraps
   for halves of 33 or more primes, and the n % L == 1 check then discards
   the genuine solutions without warning. Here no indices are stored: the
   stored-half subsets of the (rare) hits are recovered by a second pass.
3. Arithmetic overflow. The original computes (r * p) % L in uint64, which
   overflows once r * p >= 2^64 (L about 2^61 for L7). Here every residue is
   kept as a CRT pair (mod L1, mod L2) with L = L1 L2, L1, L2 < 2^32, so
   every product fits in 64 bits. The matching key is c1 + L1 * c2 < L,
   which needs L < 2^64. Larger L are skipped with a reason.

Every hit is rechecked in exact integer arithmetic: n % L == 1, and
n % (q^2 - 1) == 1 for every prime q | n (so n is a rigid CN2).

Generators are evaluated one at a time. Each evaluation already uses every
core through numba's prange, so a worker pool would only multiply memory.

Usage (edit prime_ranges in main(), as in cn2_systematic_fast.py):
    python3 cn2_systematic_fast_chunked.py [--stored-bits 32] [--max-P 70]
                                           [--threads N]
"""

import argparse
import math
import sys
import time
from typing import List, Dict, Tuple

import numpy as np
from numba import njit, prange, set_num_threads, get_num_threads

from cn2_systematic_fast import (
    divisors, calculate_P, calculate_fecundity, generate_candidates,
    factorint, GeneratorResult, print_results_table, print_latex_table,
)

LOW_BITS = 24          # streamed chunk size 2^24; low tables use 2 x 128 MiB
BUCKET_BITS = 26       # bucket index: 2^26 int64 = 512 MiB


# ============================================================================
# CRT split and residue tables
# ============================================================================

def crt_split(L_factors: Dict[int, int]) -> Tuple[int, int]:
    """Split L into coprime L1, L2 with both < 2^32, largest prime powers first.
    Returns (L1, L2) or raises ValueError."""
    powers = sorted((p ** e for p, e in L_factors.items()), reverse=True)
    L1 = L2 = 1
    for q in powers:
        if L1 <= L2 and L1 * q < 2 ** 32:
            L1 *= q
        elif L2 * q < 2 ** 32:
            L2 *= q
        elif L1 * q < 2 ** 32:
            L1 *= q
        else:
            raise ValueError("cannot split L into two factors below 2^32")
    return L1, L2


def subset_table(res1: List[int], res2: List[int], L1: int, L2: int):
    """Residues (mod L1, mod L2) of all 2^k subset products of the given
    residues. Entry j corresponds to bitmask j (bit i = element i)."""
    t1 = np.ones(1, dtype=np.uint64)
    t2 = np.ones(1, dtype=np.uint64)
    for a, b in zip(res1, res2):
        t1 = np.concatenate([t1, (t1 * np.uint64(a)) % np.uint64(L1)])
        t2 = np.concatenate([t2, (t2 * np.uint64(b)) % np.uint64(L2)])
    return t1, t2


def high_residue(h: int, res1: List[int], res2: List[int], L1: int, L2: int):
    a = b = 1
    i = 0
    while h:
        if h & 1:
            a = a * res1[i] % L1
            b = b * res2[i] % L2
        h >>= 1
        i += 1
    return a, b


# ============================================================================
# Numba kernels
# ============================================================================

@njit(parallel=True, cache=True)
def fill_keys(out, base, a, b, low1, low2, L1, L2):
    """out[base + j] = key of (a, b) * (low1[j], low2[j])."""
    n = low1.shape[0]
    for j in prange(n):
        c1 = (a * low1[j]) % L1
        c2 = (b * low2[j]) % L2
        out[base + j] = c1 + L1 * c2


@njit(parallel=True, cache=True)
def build_buckets(sorted_keys, shift, nb):
    """starts[k] = first index whose key >> shift >= k, for k = 0..nb."""
    starts = np.empty(nb + 1, dtype=np.int64)
    n = sorted_keys.shape[0]
    for k in prange(nb + 1):
        lo = np.int64(0)
        hi = np.int64(n)
        target = np.uint64(k) << np.uint64(shift)
        if k == nb:
            starts[k] = n
            continue
        while lo < hi:
            mid = (lo + hi) // 2
            if sorted_keys[mid] < target:
                lo = mid + 1
            else:
                hi = mid
        starts[k] = lo
    return starts


@njit(parallel=True, cache=True)
def stream_chunk(counts, a, b, low1, low2, L1, L2, sorted_keys, starts, shift):
    """For each j: key of (a, b) * (low1[j], low2[j]); counts[j] = number of
    stored keys equal to it."""
    n = low1.shape[0]
    for j in prange(n):
        c1 = (a * low1[j]) % L1
        c2 = (b * low2[j]) % L2
        key = c1 + L1 * c2
        bk = key >> np.uint64(shift)
        lo = starts[bk]
        hi = starts[bk + 1]
        while lo < hi:
            mid = (lo + hi) // 2
            if sorted_keys[mid] < key:
                lo = mid + 1
            else:
                hi = mid
        c = 0
        m = sorted_keys.shape[0]
        while lo < m and sorted_keys[lo] == key:
            c += 1
            lo += 1
        counts[j] = c


@njit(parallel=True, cache=True)
def mark_targets(flags, a, b, low1, low2, L1, L2, targets):
    """flags[j] = 1 if the key of (a, b) * low[j] is in the sorted targets."""
    n = low1.shape[0]
    t = targets.shape[0]
    for j in prange(n):
        c1 = (a * low1[j]) % L1
        c2 = (b * low2[j]) % L2
        key = c1 + L1 * c2
        lo = 0
        hi = t
        while lo < hi:
            mid = (lo + hi) // 2
            if targets[mid] < key:
                lo = mid + 1
            else:
                hi = mid
        flags[j] = 1 if (lo < t and targets[lo] == key) else 0


# ============================================================================
# MITM search
# ============================================================================

def mitm_search_chunked(L: int, L_factors: Dict[int, int], P: List[int],
                        stored_bits: int = 32, verbose: bool = True,
                        check_rigid: bool = True):
    """All non-empty S of P with prod(S) == 1 (mod L). Returns (hits, seconds)."""
    t0 = time.time()
    L1, L2 = crt_split(L_factors)
    assert L1 * L2 == L and math.gcd(L1, L2) == 1

    n = len(P)
    h1 = min(stored_bits, n // 2)          # stored half A = P[:h1]
    A, B = P[:h1], P[h1:]
    lo_a = min(LOW_BITS, len(A))
    lo_b = min(LOW_BITS, len(B))

    # Stored half: keys of the products.
    A1 = [p % L1 for p in A]
    A2 = [p % L2 for p in A]
    lowA1, lowA2 = subset_table(A1[:lo_a], A2[:lo_a], L1, L2)
    hiA1, hiA2 = A1[lo_a:], A2[lo_a:]
    nA = 1 << len(A)
    if verbose:
        print(f"    |A|={len(A)} stored ({nA * 8 / 2**30:.1f} GiB), "
              f"|B|={len(B)} streamed (2^{len(B)} lookups), "
              f"L1={L1}, L2={L2}", flush=True)
    keys = np.empty(nA, dtype=np.uint64)
    step = 1 << lo_a
    for h in range(1 << (len(A) - lo_a)):
        a, b = high_residue(h, hiA1, hiA2, L1, L2)
        fill_keys(keys, h * step, np.uint64(a), np.uint64(b),
                  lowA1, lowA2, np.uint64(L1), np.uint64(L2))
    keys.sort()
    shift = max(0, L.bit_length() - BUCKET_BITS)
    nb = (L - 1 >> shift) + 1
    starts = build_buckets(keys, shift, nb)
    if verbose:
        print(f"    stored half sorted in {time.time() - t0:.1f}s", flush=True)

    # Streamed half: inverse products, so a match means prod(S) == 1.
    B1 = [pow(p, -1, L1) for p in B]
    B2 = [pow(p, -1, L2) for p in B]
    lowB1, lowB2 = subset_table(B1[:lo_b], B2[:lo_b], L1, L2)
    hiB1, hiB2 = B1[lo_b:], B2[lo_b:]
    counts = np.empty(1 << lo_b, dtype=np.int64)
    raw = []                                   # (bitmask of B, key)
    nhigh = 1 << (len(B) - lo_b)
    t1 = time.time()
    for h in range(nhigh):
        a, b = high_residue(h, hiB1, hiB2, L1, L2)
        stream_chunk(counts, np.uint64(a), np.uint64(b), lowB1, lowB2,
                     np.uint64(L1), np.uint64(L2), keys, starts, shift)
        for j in np.nonzero(counts)[0]:
            j = int(j)
            c1 = a * int(lowB1[j]) % L1
            c2 = b * int(lowB2[j]) % L2
            raw.append(((h << lo_b) | j, c1 + L1 * c2))
        if verbose and nhigh >= 64 and (h + 1) % (nhigh // 64) == 0:
            el = time.time() - t1
            print(f"    streamed {h + 1}/{nhigh} chunks, {len(raw)} raw matches, "
                  f"{el:.0f}s, about {el * (nhigh - h - 1) / (h + 1):.0f}s left",
                  flush=True)
    del keys, starts

    # Recover the A-subsets for the matched keys.
    targets = np.array(sorted({k for _, k in raw}), dtype=np.uint64)
    byKey: Dict[int, List[int]] = {}
    if len(targets):
        flags = np.empty(step, dtype=np.uint8)
        for h in range(1 << (len(A) - lo_a)):
            a, b = high_residue(h, hiA1, hiA2, L1, L2)
            mark_targets(flags, np.uint64(a), np.uint64(b), lowA1, lowA2,
                         np.uint64(L1), np.uint64(L2), targets)
            for j in np.nonzero(flags)[0]:
                j = int(j)
                c1 = a * int(lowA1[j]) % L1
                c2 = b * int(lowA2[j]) % L2
                byKey.setdefault(c1 + L1 * c2, []).append((h << lo_a) | j)

    # Combine and verify exactly.
    hits = set()
    for bmB, key in raw:
        for bmA in byKey.get(key, []):
            if bmA == 0 and bmB == 0:
                continue
            S = [A[i] for i in range(len(A)) if bmA >> i & 1] + \
                [B[i] for i in range(len(B)) if bmB >> i & 1]
            nval = math.prod(S)
            if nval % L != 1:
                raise RuntimeError(f"match failed n % L == 1 for S={S}")
            if check_rigid and any(nval % (q * q - 1) != 1 for q in S):
                raise RuntimeError(f"match is not a rigid CN2: S={S}")
            hits.add(nval)
    return sorted(hits), time.time() - t0


# ============================================================================
# Generator evaluation and search
# ============================================================================

def evaluate_generator(L, L_factors, max_P_size, min_fecundity, stored_bits):
    P = calculate_P(L, L_factors, divisors(L_factors))
    k = len(P)
    fec = calculate_fecundity(L, k)

    def skipped(reason):
        return GeneratorResult(L=L, L_factors=L_factors, P_size=k, fecundity=fec,
                               cn2_count=None, min_factors=None, max_factors=None,
                               min_log10=None, max_log10=None, elapsed=None,
                               skipped=True, skip_reason=reason)
    if k > max_P_size:
        return skipped(f"|P|={k} > {max_P_size}")
    if fec < min_fecundity:
        return skipped(f"fecundity={fec:.2f} < {min_fecundity}")
    if k < 3:
        return skipped(f"|P|={k} too small")
    if L >= 2 ** 64:
        return skipped("L >= 2^64")
    try:
        crt_split(L_factors)
    except ValueError as e:
        return skipped(str(e))

    hits, el = mitm_search_chunked(L, L_factors, P, stored_bits)
    if not hits:
        return GeneratorResult(L=L, L_factors=L_factors, P_size=k, fecundity=fec,
                               cn2_count=0, min_factors=None, max_factors=None,
                               min_log10=None, max_log10=None, elapsed=el)
    nf = [len([q for q in P if h % q == 0]) for h in hits]
    return GeneratorResult(L=L, L_factors=L_factors, P_size=k, fecundity=fec,
                           cn2_count=len(hits), min_factors=min(nf),
                           max_factors=max(nf), min_log10=math.log10(hits[0]),
                           max_log10=math.log10(hits[-1]), elapsed=el), hits


def run_search(prime_ranges, outliers, max_P_size, min_fecundity,
               max_generators, stored_bits, hits_file):
    print("=" * 70)
    print("CN2 SYSTEMATIC GENERATOR SEARCH (chunked, 64-bit safe)")
    print("=" * 70)
    print(f"numba threads: {get_num_threads()}")
    cands = generate_candidates(prime_ranges, outliers)
    print(f"  {len(cands)} candidates generated")
    if max_generators:
        cands = cands[:max_generators]
    cands = [c for c in cands if c[3] <= max_P_size and c[2] >= min_fecundity]
    print(f"  {len(cands)} after filtering (|P|<={max_P_size}, fec>={min_fecundity})")

    results = []
    t0 = time.time()
    with open(hits_file, "a") as fh:
        for i, (L, Lf, fec, k) in enumerate(cands):
            print(f"  [{i + 1}/{len(cands)}] L={L}, |P|={k}, fec={fec:.2f}", flush=True)
            out = evaluate_generator(L, Lf, max_P_size, min_fecundity, stored_bits)
            if isinstance(out, tuple):
                r, hits = out
                for h in hits:
                    fh.write(f"{L}\t{h}\n")
                fh.flush()
            else:
                r = out
            if r.skipped:
                print(f"      skipped: {r.skip_reason}")
            else:
                print(f"      {r.cn2_count} CN2s in {r.elapsed:.1f}s", flush=True)
            results.append(r)

    print(f"\nTotal time: {time.time() - t0:.1f}s")
    results = [r for r in results if not r.skipped and r.cn2_count]
    results.sort(key=lambda r: -r.cn2_count)
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stored-bits", type=int, default=32,
                    help="stored half has at most this many primes "
                         "(memory 8 * 2^bits bytes; 32 -> 32 GiB, 33 -> 64 GiB)")
    ap.add_argument("--max-P", type=int, default=70)
    ap.add_argument("--min-fec", type=float, default=-20.0)
    ap.add_argument("--max-generators", type=int, default=100)
    ap.add_argument("--threads", type=int, default=0, help="numba threads (0 = all)")
    ap.add_argument("--hits-file", default="cn2_systematic_chunked_hits.tsv")
    args = ap.parse_args()
    if args.threads:
        set_num_threads(args.threads)

    # Default: the L7 generator of lean_support/Five_Seven_CN2.md (|P| = 68).
    prime_ranges = [
        (2, [8]), (3, [5]), (5, [2]), (7, [0]),
        (11, [1]), (13, [1]), (17, [1]), (19, [1]), (23, [1]),
        (29, [1]), (31, [1]), (43, [1]), (47, [1]),
    ]
    outliers = []
    print("Prime ranges:")
    for p, r in prime_ranges:
        print(f"  {p}: {r}")
    print(f"Outliers: {outliers}")

    results = run_search(prime_ranges, outliers, args.max_P, args.min_fec,
                         args.max_generators, args.stored_bits, args.hits_file)
    print_results_table(results)
    print_latex_table(results)


if __name__ == "__main__":
    main()
