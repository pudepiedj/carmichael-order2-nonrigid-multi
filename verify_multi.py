#!/usr/bin/env python3
"""verify_multi.py -- check every number in data/cn2_nonrigid_multi.tsv from scratch, using only the standard library.

    python3 verify_multi.py [data/cn2_nonrigid_multi.tsv] [--afile data/a175531_multi.txt] [--searches data/searches/*.jsonl]
        [--ktree data/cn2_ktree_multi.tsv --ktree-searches data/searches/cn2_ktree_runs.jsonl]
        [--div57 data/cn2_divisible_5_7.tsv]

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
With --ktree: a second dataset in the same format (the Wagner k-tree finds), checked row by row in the same way.
With --ktree-searches: the k-tree run records (L0, p0s, count, values); every value must be in the k-tree dataset with
type_p == p0s and n == 1 (mod L0), times_found must equal the number of runs reporting n, and every k-tree row must
come from some run.
With --div57: the numbers divisible by 5 or 7 (columns log10 k divisor generator method n primes). For every row
the primes are checked as above, and n == 1 (mod p^2 - 1) at EVERY prime p | n (all are rigid), the divisor is 5
or 7 and divides n, and the rows are strictly increasing in n within each divisor.
With --afile: every number in the a-file is in one of the datasets given, with the same prime factors and the same
primes p (separated by commas, or by '.' in the version-1.0 a-file; '-' means none).
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


def check_div_row(f):
    """-> (n, divisor, None) if the row is good, else (n or None, None, reason)."""
    try:
        if len(f) != 7:
            return None, None, f'{len(f)} fields'
        lg, k, dv, _gen, _method, n, ps = f
        n, dv = int(n), int(dv)
        primes = [int(x) for x in ps.split()]
        if len(set(primes)) != len(primes) or primes != sorted(primes):
            return n, None, 'primes repeated or not sorted'
        if not all(p < 3 * 10**24 and is_prime(p) for p in primes):
            return n, None, 'a listed factor is not prime'
        if int(k) != len(primes) or len(primes) < 2:
            return n, None, 'k is not the number of primes'
        if math.prod(primes) != n:
            return n, None, 'product of primes is not n'
        if f'{math.log10(n):.3f}' != lg:
            return n, None, f'log10 {lg} != {math.log10(n):.3f}'
        if dv not in (5, 7) or dv not in primes:
            return n, None, f'divisor {dv} is not 5 or 7, or does not divide n'
        for p in primes:
            if n % (p * p - 1) != 1:
                return n, None, f'not rigid at p = {p}'
        return n, dv, None
    except Exception as e:
        return None, None, f'unparseable: {e}'


def read_tsv(path, header):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        h = fh.readline().rstrip('\n').split('\t')
        rows = [line.rstrip('\n').split('\t') for line in fh if line.strip()]
    if h != header:
        sys.exit(f'{path}: unexpected header {h}')
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('data', nargs='?', default='data/cn2_nonrigid_multi.tsv')
    ap.add_argument('--afile')
    ap.add_argument('--searches', nargs='+', default=[])
    ap.add_argument('--ktree')
    ap.add_argument('--ktree-searches', nargs='+', default=[])
    ap.add_argument('--div57')
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

    kfail = 0
    factors = {int(f[5]): sorted(int(x) for x in f[6].split()) for f in rows}
    if a.ktree:
        krows = read_tsv(a.ktree, HEADER)
        kres = [check_row(f) for f in krows]
        prev = 0
        for i, (n, t, c, why) in enumerate(kres):
            if why or n <= prev or n in good:
                kfail += 1
                print(f'K-TREE FAIL line {i + 2}: n = {str(n)[:30]}: {why or "out of order or already in the main dataset"}')
            prev = n or prev
        kgood = {n: t for n, t, _, why in kres if not why}
        ktimes = {n: c for n, _, c, why in kres if not why}
        for f in krows:
            factors[int(f[5])] = sorted(int(x) for x in f[6].split())
        kc = Counter(len(t) for t in kgood.values())
        print(f'k-tree: {len(krows)} rows, {len(kgood)} verified; by number of primes p: '
              + ', '.join(f'{k}: {v}' for k, v in sorted(kc.items())))
        if a.ktree_searches:
            reported = defaultdict(int)
            runs = 0
            for path in a.ktree_searches:
                for line in open(path):
                    if not line.strip():
                        continue
                    r = json.loads(line)
                    runs += 1
                    L0, T = int(r['L0']), tuple(sorted(int(x) for x in r['p0s']))
                    vals = [int(v) for v in r['values']]
                    if len(vals) != int(r['count']):
                        kfail += 1
                        print(f'K-TREE RUN {T} run {r.get("run")}: count {r["count"]} but {len(vals)} values')
                    for n in set(vals):
                        reported[n] += 1
                        if kgood.get(n) != T or n % L0 != 1:
                            kfail += 1
                            print(f'K-TREE RUN {T}: value {str(n)[:30]}... not in the k-tree dataset with these primes p, or n != 1 mod L0')
            for n in kgood:
                if reported.get(n, 0) != ktimes[n]:
                    kfail += 1
                    print(f'K-TREE TIMES: {str(n)[:30]}... times_found {ktimes[n]} but {reported.get(n, 0)} runs report it')
            print(f'k-tree runs: {runs} records, {sum(reported.values())} reports, every k-tree row accounted for: {all(n in reported for n in kgood)}')
        good.update(kgood)

    dfail = 0
    if a.div57:
        drows = read_tsv(a.div57, ['log10', 'k', 'divisor', 'generator', 'method', 'n', 'primes'])
        dres = [check_div_row(f) for f in drows]
        last = {}
        for i, (n, dv, why) in enumerate(dres):
            if not why and (n <= last.get(dv, 0) or n in good):
                why = 'out of order within its divisor, or already listed'
            if why:
                dfail += 1
                print(f'DIV57 FAIL line {i + 2}: n = {str(n)[:30]}: {why}')
            else:
                last[dv] = n
                good[n] = ()
                factors[n] = sorted(int(x) for x in drows[i][6].split())
        dc = Counter(dv for _, dv, why in dres if not why)
        print(f'divisible by 5 or 7: {len(drows)} rows, {sum(dc.values())} verified (5: {dc[5]}, 7: {dc[7]}), '
              f'smallest with 5: {min((len(str(n)) for n, dv, w in dres if dv == 5), default=0)} digits, '
              f'with 7: {min((len(str(n)) for n, dv, w in dres if dv == 7), default=0)} digits')

    afail = 0
    if a.afile:
        seen = 0
        for line in open(a.afile):
            if not line.strip() or line.startswith('#'):
                continue
            t = line.split()
            n = int(t[0])
            tset = () if t[2] == '-' else tuple(sorted(int(x) for x in t[2].replace(',', '.').split('.')))   # v1.0 a-file uses '.', v1.1 ','
            seen += 1
            if good.get(n) != tset or (len(t) > 3 and sorted(int(x) for x in t[3:]) != factors.get(n)):
                afail += 1
                print(f'A-FILE: {t[0][:30]}... not in the datasets with primes p = {t[2]} and the same factors')
        print(f'a-file: {seen} numbers, {seen - afail} found in the datasets with the same primes p and factors')

    ok = not bad and not order and not afail and not sfail and not kfail and not dfail
    print('ALL ROWS VERIFIED' if ok else 'VERIFICATION FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
