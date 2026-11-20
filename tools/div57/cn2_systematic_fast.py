#!/usr/bin/env python3
"""
CN2 Systematic Generator Search - Parallelized

Explores L = 2^k1 × 3^k2 × ... × p^kn systematically,
using fast MITM to evaluate each generator.

Features:
- Parallel evaluation of generators
- Best-first search by estimated fecundity
- LaTeX table output
- Flexible prime ranges and outliers
"""

import numpy as np
from numba import njit, prange
import math
import time
import itertools
from sympy import isprime, totient
from functools import reduce
from operator import mul
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
import sys

# ============================================================================
# CORE MATH (same as cn2_fast.py)
# ============================================================================

def factorint(n):
    d, out = 2, {}
    while d*d <= n:
        while n % d == 0:
            out[d] = out.get(d,0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1: out[n] = out.get(n,0) + 1
    return out

def divisors(fact):
    divs = [1]
    for p, e in fact.items():
        divs = [d * p**k for d in divs for k in range(e+1)]
    return divs

def calculate_P(L, factL, divs):
    P = []
    for d in divs:
        s = d + 1
        root = math.isqrt(s)
        if root*root == s and root > 1 and isprime(root) and root not in factL:
            P.append(root)
    return sorted(P)

def calculate_fecundity(L, P_size):
    """Fecundity = |P| - log2(φ(L))"""
    phi_L = totient(L)
    return P_size - math.log2(phi_L)

# ============================================================================
# NUMBA-ACCELERATED MITM (from cn2_fast.py)
# ============================================================================

@njit(cache=True)
def mod_inverse(a, m):
    if a == 0:
        return np.uint64(0)
    m0, x0, x1 = m, np.int64(0), np.int64(1)
    a = a % m
    while a > 1:
        q = a // m
        m, a = a % m, m
        x0, x1 = x1 - q * x0, x0
    if x1 < 0:
        x1 += m0
    return np.uint64(x1)

@njit(parallel=True, cache=True)
def compute_residues_parallel(primes, L, n_bits):
    n_subsets = 1 << n_bits
    residues = np.empty(n_subsets, dtype=np.uint64)
    for bitmask in prange(n_subsets):
        r = np.uint64(1)
        for i in range(n_bits):
            if bitmask & (1 << i):
                r = (r * np.uint64(primes[i])) % np.uint64(L)
        residues[bitmask] = r
    return residues

@njit(parallel=True, cache=True)
def find_matches_parallel(residues2, sorted_residues1, sort_indices1, L):
    n2 = len(residues2)
    n1 = len(sorted_residues1)
    
    match_counts = np.zeros(n2, dtype=np.int64)
    
    for i in prange(n2):
        r2 = residues2[i]
        if r2 == 0:
            continue
        target = mod_inverse(r2, L)
        if target == 0:
            continue
        lo, hi = np.int64(0), np.int64(n1)
        while lo < hi:
            mid = (lo + hi) // 2
            if sorted_residues1[mid] < target:
                lo = mid + 1
            else:
                hi = mid
        count = 0
        while lo < n1 and sorted_residues1[lo] == target:
            count += 1
            lo += 1
        match_counts[i] = count
    
    total_matches = np.sum(match_counts)
    offsets = np.empty(n2, dtype=np.int64)
    offset = np.int64(0)
    for i in range(n2):
        offsets[i] = offset
        offset += match_counts[i]
    
    matches = np.empty((total_matches, 2), dtype=np.uint32)
    
    for i in prange(n2):
        r2 = residues2[i]
        if r2 == 0:
            continue
        target = mod_inverse(r2, L)
        if target == 0:
            continue
        lo, hi = np.int64(0), np.int64(n1)
        while lo < hi:
            mid = (lo + hi) // 2
            if sorted_residues1[mid] < target:
                lo = mid + 1
            else:
                hi = mid
        write_idx = offsets[i]
        while lo < n1 and sorted_residues1[lo] == target:
            matches[write_idx, 0] = i
            matches[write_idx, 1] = sort_indices1[lo]
            write_idx += 1
            lo += 1
    
    return matches

def bitmask_to_product(bitmask, primes):
    prod = 1
    for i, p in enumerate(primes):
        if bitmask & (1 << i):
            prod *= int(p)
    return prod

# ============================================================================
# FAST MITM SEARCH (streamlined for batch use)
# ============================================================================

def mitm_search(L: int, P: List[int], P_prod: int, verbose: bool = False) -> Tuple[List[int], float]:
    """
    Run MITM search for CN2s.
    Returns (list of CN2s, elapsed time in seconds)
    """
    n = len(P)
    half = n // 2
    bag1 = np.array(P[:half], dtype=np.int64)
    bag2 = np.array(P[half:], dtype=np.int64)
    
    t0 = time.time()
    
    # Compute residues
    residues1 = compute_residues_parallel(bag1, L, len(bag1))
    residues2 = compute_residues_parallel(bag2, L, len(bag2))
    
    # Sort bag1
    sort_indices1 = np.argsort(residues1).astype(np.uint32)
    sorted_residues1 = residues1[sort_indices1]
    
    # Find matches
    matches = find_matches_parallel(residues2, sorted_residues1, sort_indices1, L)
    
    # Verify matches
    hits = []
    bag1_list = [int(x) for x in bag1]
    bag2_list = [int(x) for x in bag2]
    
    for i in range(len(matches)):
        bm2 = int(matches[i, 0])
        bm1 = int(matches[i, 1])
        
        if bm1 == 0 and bm2 == 0:
            continue
        
        p1 = bitmask_to_product(bm1, bag1_list)
        p2 = bitmask_to_product(bm2, bag2_list)
        n_val = p1 * p2
        
        if n_val == 1:
            continue
        
        if n_val % L != 1:
            continue
        
        if P_prod % n_val == 0:
            hits.append(n_val)
    
    hits = sorted(set(hits))
    elapsed = time.time() - t0
    
    return hits, elapsed

# ============================================================================
# GENERATOR EVALUATION
# ============================================================================

@dataclass
class GeneratorResult:
    """Results from evaluating a generator"""
    L: int
    L_factors: Dict[int, int]
    P_size: int
    fecundity: float
    cn2_count: Optional[int]
    min_factors: Optional[int]
    max_factors: Optional[int]
    min_log10: Optional[float]
    max_log10: Optional[float]
    elapsed: Optional[float]
    skipped: bool = False
    skip_reason: str = ""

def evaluate_generator(L: int, L_factors: Dict[int, int], 
                       max_P_size: int = 66,
                       min_fecundity: float = -20.0) -> GeneratorResult:
    """
    Evaluate a single generator L.
    Returns GeneratorResult with all metrics.
    """
    # Calculate P
    divs = divisors(L_factors)
    P = calculate_P(L, L_factors, divs)
    P_size = len(P)
    
    # Calculate fecundity estimate
    fecundity = calculate_fecundity(L, P_size)
    
    # Check skip conditions
    if P_size > max_P_size:
        return GeneratorResult(
            L=L, L_factors=L_factors, P_size=P_size, fecundity=fecundity,
            cn2_count=None, min_factors=None, max_factors=None,
            min_log10=None, max_log10=None, elapsed=None,
            skipped=True, skip_reason=f"|P|={P_size} > {max_P_size}"
        )
    
    if fecundity < min_fecundity:
        return GeneratorResult(
            L=L, L_factors=L_factors, P_size=P_size, fecundity=fecundity,
            cn2_count=None, min_factors=None, max_factors=None,
            min_log10=None, max_log10=None, elapsed=None,
            skipped=True, skip_reason=f"fecundity={fecundity:.2f} < {min_fecundity}"
        )
    
    if P_size < 3:
        return GeneratorResult(
            L=L, L_factors=L_factors, P_size=P_size, fecundity=fecundity,
            cn2_count=0, min_factors=None, max_factors=None,
            min_log10=None, max_log10=None, elapsed=0.0,
            skipped=True, skip_reason=f"|P|={P_size} too small"
        )
    
    # Compute P_prod
    P_prod = reduce(mul, P, 1)
    
    # Run MITM
    hits, elapsed = mitm_search(L, P, P_prod)
    
    if len(hits) == 0:
        return GeneratorResult(
            L=L, L_factors=L_factors, P_size=P_size, fecundity=fecundity,
            cn2_count=0, min_factors=None, max_factors=None,
            min_log10=None, max_log10=None, elapsed=elapsed
        )
    
    # Analyze hits
    factor_counts = [len(factorint(h)) for h in hits]
    
    return GeneratorResult(
        L=L, L_factors=L_factors, P_size=P_size, fecundity=fecundity,
        cn2_count=len(hits),
        min_factors=min(factor_counts),
        max_factors=max(factor_counts),
        min_log10=math.log10(hits[0]),
        max_log10=math.log10(hits[-1]),
        elapsed=elapsed
    )

# ============================================================================
# GENERATOR ENUMERATION
# ============================================================================

def build_L_from_exponents(primes: List[int], exponents: Tuple[int, ...]) -> Tuple[int, Dict[int, int]]:
    """Build L and L_factors from primes and exponents"""
    L = 1
    L_factors = {}
    for prime, exp in zip(primes, exponents):
        if exp > 0:
            L *= prime ** exp
            L_factors[prime] = exp
    return L, L_factors

def generate_candidates(prime_ranges: List[Tuple[int, List[int]]], 
                        outliers: List[int] = None) -> List[Tuple[int, Dict[int, int], float]]:
    """
    Generate all candidate L values with estimated fecundity.
    Returns list of (L, L_factors, fecundity_estimate) sorted by fecundity descending.
    """
    if outliers is None:
        outliers = []
    
    # Add outliers to prime_ranges with exponent [0, 1]
    full_ranges = list(prime_ranges)
    for o in outliers:
        full_ranges.append((o, [0, 1]))
    
    primes = [p for p, _ in full_ranges]
    ranges = [r for _, r in full_ranges]
    
    candidates = []
    
    for exponents in itertools.product(*ranges):
        L, L_factors = build_L_from_exponents(primes, exponents)
        
        if L == 1:
            continue
        
        # Quick P size estimate (actually calculate it)
        divs = divisors(L_factors)
        P = calculate_P(L, L_factors, divs)
        P_size = len(P)
        
        if P_size < 3:
            continue
        
        fecundity = calculate_fecundity(L, P_size)
        candidates.append((L, L_factors, fecundity, P_size))
    
    # Sort by fecundity (descending) - best first
    candidates.sort(key=lambda x: -x[2])
    
    return candidates

# ============================================================================
# PARALLEL SEARCH
# ============================================================================

def run_parallel_search(prime_ranges: List[Tuple[int, List[int]]],
                        outliers: List[int] = None,
                        max_P_size: int = 66,
                        min_fecundity: float = -20.0,
                        max_generators: int = None,
                        n_workers: int = None,
                        verbose: bool = True) -> List[GeneratorResult]:
    """
    Run parallel search over generator space.
    
    Args:
        prime_ranges: List of (prime, [exponent options])
        outliers: Additional primes to include with exponent [0,1]
        max_P_size: Skip generators with |P| > this
        min_fecundity: Skip generators with fecundity < this
        max_generators: Only evaluate top N by estimated fecundity
        n_workers: Number of parallel workers (default: CPU count)
        verbose: Print progress
    
    Returns:
        List of GeneratorResult sorted by cn2_count descending
    """
    if outliers is None:
        outliers = []
    
    if n_workers is None:
        n_workers = min(8, len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else 8)
    
    # Generate candidates
    if verbose:
        print("="*70)
        print("CN2 SYSTEMATIC GENERATOR SEARCH")
        print("="*70)
        print(f"Generating candidate L values...")
    
    candidates = generate_candidates(prime_ranges, outliers)
    
    if verbose:
        print(f"  {len(candidates)} candidates generated")
        print(f"  Top 5 by estimated fecundity:")
        for i, (L, factors, fec, P_size) in enumerate(candidates[:5]):
            print(f"    {i+1}. L={L}, |P|={P_size}, fecundity≈{fec:.2f}")
    
    # Limit if requested
    if max_generators and len(candidates) > max_generators:
        candidates = candidates[:max_generators]
        if verbose:
            print(f"  Limiting to top {max_generators} candidates")
    
    # Filter by constraints
    filtered = [(L, factors, fec, P_size) for L, factors, fec, P_size in candidates 
                if P_size <= max_P_size and fec >= min_fecundity]
    
    if verbose:
        print(f"  {len(filtered)} candidates after filtering (|P|≤{max_P_size}, fec≥{min_fecundity})")
    
    if len(filtered) == 0:
        print("No candidates to evaluate!")
        return []
    
    # Evaluate in parallel
    if verbose:
        print(f"\nEvaluating with {n_workers} workers...")
    
    results = []
    t0 = time.time()
    
    # For small numbers, just run sequentially to avoid overhead
    if len(filtered) <= 4 or n_workers == 1:
        for i, (L, L_factors, _, _) in enumerate(filtered):
            result = evaluate_generator(L, L_factors, max_P_size, min_fecundity)
            results.append(result)
            if verbose and not result.skipped:
                print(f"  [{i+1}/{len(filtered)}] L={L}: {result.cn2_count} CN2s in {result.elapsed:.2f}s")
    else:
        # Parallel evaluation
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(evaluate_generator, L, L_factors, max_P_size, min_fecundity): (L, L_factors)
                for L, L_factors, _, _ in filtered
            }
            
            completed = 0
            for future in as_completed(futures):
                completed += 1
                result = future.result()
                results.append(result)
                
                if verbose and not result.skipped:
                    print(f"  [{completed}/{len(filtered)}] L={result.L}: "
                          f"{result.cn2_count} CN2s, |P|={result.P_size}, "
                          f"fec={result.fecundity:.2f}, {result.elapsed:.2f}s")
    
    total_time = time.time() - t0
    
    # Sort by cn2_count descending
    results = [r for r in results if not r.skipped and r.cn2_count and r.cn2_count > 0]
    results.sort(key=lambda r: -r.cn2_count)
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"SEARCH COMPLETE")
        print(f"{'='*70}")
        print(f"  Total time: {total_time:.1f}s")
        print(f"  Generators with CN2s: {len(results)}")
        if results:
            print(f"  Most productive: L={results[0].L} with {results[0].cn2_count} CN2s")
    
    return results

# ============================================================================
# OUTPUT FORMATTING
# ============================================================================

def format_L_factors(L_factors: Dict[int, int]) -> str:
    """Format L factors as 2^7.3^3.5^2..."""
    parts = []
    for p in sorted(L_factors.keys()):
        e = L_factors[p]
        if e == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{e}")
    return '.'.join(parts)

def print_results_table(results: List[GeneratorResult], top_n: int = 20):
    """Print results as a formatted table"""
    print(f"\n{'='*100}")
    print(f"TOP {min(top_n, len(results))} GENERATORS BY CN2 COUNT")
    print(f"{'='*100}")
    print(f"{'Rank':>4} {'L':>20} {'|P|':>4} {'Fec':>6} {'CN2s':>8} {'MinF':>4} {'MaxF':>4} {'log₁₀ range':>14} {'Time':>8}")
    print("-"*100)
    
    for i, r in enumerate(results[:top_n]):
        log_range = f"{r.min_log10:.1f}-{r.max_log10:.1f}" if r.min_log10 else "N/A"
        print(f"{i+1:>4} {r.L:>20} {r.P_size:>4} {r.fecundity:>6.2f} {r.cn2_count:>8} "
              f"{r.min_factors or 0:>4} {r.max_factors or 0:>4} {log_range:>14} {r.elapsed:>7.2f}s")

def print_latex_table(results: List[GeneratorResult], top_n: int = 20):
    """Print results as LaTeX table rows"""
    print(f"\n% LaTeX table rows (top {min(top_n, len(results))})")
    print(r"% L & $L_{factors}$ & CN2s & |P| & MinF & MaxF & log_{min} & log_{max} & time \\")
    
    for r in results[:top_n]:
        factors_str = format_L_factors(r.L_factors)
        print(f"L & ${r.L}$ & ${factors_str}$ & ${r.cn2_count}$ & ${r.P_size}$ & "
              f"{r.min_factors} & {r.max_factors} & {r.min_log10:.2f} & {r.max_log10:.2f} & "
              f"{r.elapsed:.4f}s \\\\")

# ============================================================================
# MAIN
# ============================================================================

import os

def main():
    # Default prime ranges (adjust as needed)
    prime_ranges = [
        (2,  [8]),
        (3,  [5]),
        (5,  [2]),
        (7,  [0]),
        (11, [1]),
        (13, [1]),
        (17, [1]),
        (19, [1]),
        (23, [1]),
        (29, [1]),
        (31, [1]),
        (43, [1]),
        (47, [1])
    ]
    
    # Outliers to test
    outliers = []
    
    print("Prime ranges:")
    for p, r in prime_ranges:
        print(f"  {p}: {r}")
    print(f"Outliers: {outliers}")
    
    # Run search
    results = run_parallel_search(
        prime_ranges=prime_ranges,
        outliers=outliers,
        max_P_size=70,       # Adjust based on time budget
        min_fecundity=-20.0,  # below corpus floor (~-13.2): keep deep-tail sporadics
        max_generators=100,  # Limit to top 100 by estimated fecundity
        n_workers=8,
        verbose=True
    )
    
    # Output
    print_results_table(results)
    print_latex_table(results)
    
    return results


if __name__ == "__main__":
    results = main()