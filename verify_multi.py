#!/usr/bin/env python3
"""verify_multi.py -- check every number in data/cn2_nonrigid_multi.tsv from scratch, using only the standard library.

    python3 verify_multi.py [data/cn2_nonrigid_multi.tsv] [--afile data/a175531_multi.txt] [--searches data/searches/*.jsonl]

The dataset is tab-separated (gzip-compressed if the name ends in .gz), one header line, then one row per number,
sorted by n:

    log10  k  type_p  L_R  times_found  n  primes

For every row this checks that:
  - the listed primes are distinct and each is prime (deterministic Miller-Rabin; valid below 3.3e24);
  - k is the number of primes, their product is n, and log10 is log10(n) to 3 decimals;
  - n is a Carmichael number of order 2 (Howe's definition): n == 1 or p (mod p^2 - 1) for every prime p | n;
  - the primes p with n == p (not 1) mod p^2 - 1 are exactly the listed type_p, and there are at least two;
  - L_R is lcm(q^2 - 1) over the other primes q;
  - the rows are strictly increasing in n, so no number is listed twice.
With --afile: every number in the a-file is in the dataset with the same primes p.
With --searches: the search records (one JSON object per line: L0, p0s, predicted, count, values).  Each distinct
search (L0, p0s) is counted once; every value it reports must be in the dataset with type_p == p0s and n == 1
(mod L0); times_found must equal the number of distinct searches reporting n; and every row must come from at
least one search.  It then prints found against predicted for two and for three primes.
Prints a summary and exits non-zero on any failure.
"""
import argparse, gzip, json, math, sys
from collections import Counter, defaultdict

MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)   # deterministic for n < 3.3e24
HEADER = ['log10', 'k', 'type_p', 'L_R', 'times_found', 'n', 'primes']


def is_prime(n):
    if n < 2:
        return False
    for p in MR_BASES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in MR_BASES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def check_row(f):
    """-> (n, type_p tuple, times, None) if the row is good, else (n or None, None, None, reason)."""
    try:
        if len(f) != 7:
            return None, None, None, f'{len(f)} fields'
        lg, k, tp, lr, times, n, ps = f
        n = int(n)
        primes = [int(x) for x in ps.split()]
        tset = tuple(sorted(int(x) for x in tp.split(',')))
        if len(set(primes)) != len(primes):
            return n, None, None, 'repeated prime'
        if not all(p < 3 * 10**24 and is_prime(p) for p in primes):
            return n, None, None, 'a listed factor is not prime'
        if int(k) != len(primes):
            return n, None, None, 'k is not the number of primes'
        if math.prod(primes) != n:
            return n, None, None, 'product of primes is not n'
        if f'{math.log10(n):.3f}' != lg:
            return n, None, None, f'log10 {lg} != {math.log10(n):.3f}'
        typep, L = [], 1
        for p in primes:
            r = n % (p * p - 1)
            if r == p:
                typep.append(p)
            elif r == 1:
                L = math.lcm(L, p * p - 1)
            else:
                return n, None, None, f'not of order 2 at p = {p}'
        if len(typep) < 2:
            return n, None, None, f'only {len(typep)} prime(s) with n == p (mod p^2-1)'
        if tuple(sorted(typep)) != tset:
            return n, None, None, f'primes {sorted(typep)} != listed {tp}'
        if L != int(lr):
            return n, None, None, 'L_R is not lcm(q^2-1) over the other primes'
        if int(times) < 1:
            return n, None, None, 'times_found < 1'
        return n, tset, int(times), None
    except Exception as e:                                        # malformed row
        return None, None, None, f'unparseable: {e}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('data', nargs='?', default='data/cn2_nonrigid_multi.tsv')
    ap.add_argument('--afile')
    ap.add_argument('--searches', nargs='+', default=[])
    a = ap.parse_args()
    op = gzip.open if a.data.endswith('.gz') else open
    with op(a.data, 'rt') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        rows = [line.rstrip('\n').split('\t') for line in fh if line.strip()]
    if header != HEADER:
        sys.exit(f'unexpected header {header}')
    res = [check_row(f) for f in rows]
    bad = [(i + 2, n, why) for i, (n, _, _, why) in enumerate(res) if why]
    prev, order = 0, 0
    for n, _, _, why in res:
        if n is not None and not why:
            order += n <= prev
            prev = n
    good = {n: t for n, t, _, why in res if not why}
    times = {n: c for n, _, c, why in res if not why}
    by_count = Counter(len(t) for t in good.values())
    summary = {'rows': len(rows), 'verified': len(good), 'failed': len(bad), 'order_violations': order,
               'by_number_of_primes_p': {str(k): v for k, v in sorted(by_count.items())}}
    for k in sorted(by_count):
        ns = [n for n, t in good.items() if len(t) == k]
        summary[f'smallest_with_{k}'] = str(min(ns))
        summary[f'log10_range_with_{k}'] = [round(math.log10(min(ns)), 3), round(math.log10(max(ns)), 3)]
    print(json.dumps(summary, indent=1))
    for line, n, why in bad[:20]:
        print(f'FAIL line {line}: n = {n}: {why}')

    afail = 0
    if a.afile:
        seen = 0
        for line in open(a.afile):
            if not line.strip() or line.startswith('#'):
                continue
            t = line.split()
            n, tset = int(t[0]), tuple(sorted(int(x) for x in t[2].split('.')))
            seen += 1
            if good.get(n) != tset:
                afail += 1
                print(f'A-FILE: {t[0][:30]}... not in the dataset with primes p = {t[2]}')
        print(f'a-file: {seen} numbers, {seen - afail} found in the dataset with the same primes p')

    sfail = 0
    if a.searches:
        searches = {}
        for path in a.searches:
            for line in open(path):
                if line.strip():
                    r = json.loads(line)
                    searches[(int(r['L0']), tuple(sorted(r['p0s'])))] = r     # a repeated search is counted once
        reported = defaultdict(int)
        tally = defaultdict(lambda: [0, 0.0, 0, set()])                      # searches, predicted, hits, distinct
        for (L0, T), r in searches.items():
            vals = [int(v) for v in r['values']]
            if len(vals) != r['count']:
                sfail += 1
                print(f'SEARCH {L0} {T}: count {r["count"]} but {len(vals)} values')
            t = tally[len(T)]
            t[0] += 1; t[1] += r['predicted']; t[2] += len(vals); t[3].update(vals)
            for n in set(vals):
                reported[n] += 1
                if good.get(n) != T or n % L0 != 1:
                    sfail += 1
                    print(f'SEARCH {L0} {T}: value {str(n)[:30]}... not in the dataset with these primes p, or n != 1 mod L0')
        for n in good:
            if reported.get(n, 0) != times[n]:
                sfail += 1
                print(f'TIMES: {str(n)[:30]}... times_found {times[n]} but {reported.get(n, 0)} searches report it')
        for k in sorted(tally):
            s, pred, hits, distinct = tally[k]
            print(f'{k} primes p: {s} searches, {hits} hits against {pred:.1f} predicted, {len(distinct)} distinct numbers')

    ok = not bad and not order and not afail and not sfail
    print('ALL ROWS VERIFIED' if ok else 'VERIFICATION FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
