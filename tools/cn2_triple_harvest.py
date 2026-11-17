#!/usr/bin/env python3
"""cn2_triple_harvest.py -- run the priced multi-type-p sets of cn2_multi_scout.py, best yield per hour first.

    python3 cn2_multi_scout.py --k 3 --fmin 24 --pmax 92 --top 14 --emin 0.003 --out cn2_triple_sets.json
    python3 cn2_triple_harvest.py --pmax 84 --hours 8          # resumable: repeat the same command
    python3 cn2_triple_harvest.py --report
    # two type-p primes: the pairs of cn2_twin_scout.py not yet run (records in their own file)
    python3 cn2_triple_harvest.py --sets cn2_twin_pairs.json --out cn2_twin_harvest.jsonl --skip cn2_twin_runs.jsonl \
                                  --pmax 72 --emin 0.05 --hours 14 --wait

Queue order: expected count E divided by estimated run time (~2^(|P'|/2) list entries on cn2ss4), restricted to
pools <= --pmax.  Every finished (L0, T) is one line of cn2_triple_harvest.jsonl and is skipped on a rerun; an
interrupted engine run continues from its own checkpoint (Ctrl-C once, wait for the engine to stop).
The same number can turn up from several generators, so --report counts DISTINCT numbers.

Only one harvest runs at a time (a lock file); --wait queues this one behind the harvest now running.
Before anything is run, every set is checked against the rules of note section 2 (check_rules): all members
== 1 mod 24 (rule 3), equal v2(p-1) (rule 2), no odd prime dividing p+1 of one member and p'-1 of another, and no
member dividing another's p^2-1 (rule 4).  Rule 1 (private primes) concerns the generator and is enforced by the
admissibility test in search_multi_fast.  A set that fails is a scout error and stops the harvest.
"""
import argparse, json, math, os, time

from sympy import factorint

from cn2_nonrigid_search_fast import search_multi_fast

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'cn2_triple_harvest.jsonl')


def v2(x):
    return (x & -x).bit_length() - 1


def check_rules(T):
    """None if T satisfies rules 2-4 of the note, else the reason it fails."""
    for p in T:
        if p % 24 != 1:
            return f"rule 3: {p} is not 1 mod 24"
    if len({v2(p - 1) for p in T}) != 1:
        return "rule 2: members differ in v2(p-1)"
    for p in T:
        for r in T:
            if p == r:
                continue
            g = math.gcd(p + 1, r - 1)
            if g >> v2(g) != 1:
                return f"rule 4: an odd prime divides {p}+1 and {r}-1"
            if (r * r - 1) % p == 0:
                return f"rule 4: {p} divides {r}^2-1"
    return None


def est_seconds(P, threads=32):
    half = P // 2
    return (2.0 ** half + 2.0 ** (P - half)) * (8e-9 if P <= 86 else 13e-9) / threads + (2 if P < 72 else 5 if P < 78 else 25)


def load_done(path=None):
    """one record per (L0, T): a repeated line (two harvests run at once) is counted once"""
    done = {}
    path = path or OUT
    if os.path.exists(path):
        for line in open(path):
            if line.strip():
                r = json.loads(line)
                done.setdefault((str(r['L0']), tuple(r['p0s'])), r)
    return list(done.values())


def report():
    done = load_done()
    vals = {}
    for r in done:
        for v in r['values']:
            vals.setdefault(int(v), []).append((r['L0'], r['p0s']))
    print(f"{len(done)} runs, {sum(r['seconds'] for r in done) / 3600:.2f} h, predicted {sum(r['predicted'] for r in done):.2f}, "
          f"found {sum(r['count'] for r in done)} ({len(vals)} distinct), all verified: "
          f"{all(r['verified'] == r['count'] for r in done)}")
    for n in sorted(vals):
        print(f"  10^{math.log10(n):7.2f}  T={vals[n][0][1]}  from {len(vals[n])} generator(s)  n = {n}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sets', default=os.path.join(HERE, 'cn2_triple_sets.json'))
    ap.add_argument('--pmax', type=int, default=84)
    ap.add_argument('--hours', type=float, default=8.0, help='stop starting new runs after this long')
    ap.add_argument('--max-runs', type=int, default=10 ** 9)
    ap.add_argument('--threads', type=int, default=32)
    ap.add_argument('--report', action='store_true')
    ap.add_argument('--out', help='record file (default cn2_triple_harvest.jsonl)')
    ap.add_argument('--skip', help='records of an earlier runner; its (L0, T) are not repeated')
    ap.add_argument('--emin', type=float, default=0.0)
    ap.add_argument('--dry-run', action='store_true', help='show the queue and stop')
    ap.add_argument('--wait', action='store_true', help='wait for the harvest now running to end, then start')
    a = ap.parse_args()
    global OUT
    if a.out:
        OUT = os.path.join(HERE, a.out)
    if a.report:
        return report()
    sets = [d for d in json.load(open(a.sets)) if d['P'] <= a.pmax and d['E'] >= a.emin]
    for d in sets:
        d.setdefault('T', d.get('p0s'))                       # cn2_twin_pairs.json calls it p0s
        why = check_rules(d['T'])
        if why:
            raise SystemExit(f"{d['T']} on L0 = {d['L0']} breaks {why}")
    sets.sort(key=lambda d: -d['E'] / est_seconds(d['P'], a.threads))
    done = {(str(r['L0']), tuple(r['p0s'])) for r in load_done()}
    if a.skip:
        done |= {(str(r['L0']), tuple(r['p0s'])) for r in load_done(os.path.join(HERE, a.skip))}
    todo = [d for d in sets if (d['L0'], tuple(d['T'])) not in done]
    print(f"{len(sets)} sets with pool <= {a.pmax}, all within rules 2-4; {len(sets) - len(todo)} done, {len(todo)} to go "
          f"(E {sum(d['E'] for d in todo):.1f}, ~{sum(est_seconds(d['P'], a.threads) for d in todo) / 3600:.1f} h)", flush=True)
    if a.dry_run:
        for d in todo[:10]:
            print(f"  pool {d['P']}  E = {d['E']:.3f}  est {est_seconds(d['P'], a.threads):.0f}s  T = {d['T']}  L0 = {d['L0']}")
        return
    import fcntl
    lock = open(os.path.join(HERE, 'cn2_triple_harvest.lock'), 'w')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)      # held until this process exits
    except OSError:
        if not a.wait:
            raise SystemExit("another cn2_triple_harvest.py is already running here; two would only repeat each "
                             "other's work (--wait queues this one behind it)")
        print("waiting for the harvest now running to end ...", flush=True)
        fcntl.flock(lock, fcntl.LOCK_EX)
        print(f"started {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    t0, runs = time.time(), 0
    for d in todo:
        if time.time() - t0 > a.hours * 3600 or runs >= a.max_runs:
            break
        L0 = int(d['L0'])
        print(f"\n[{runs + 1}] L0 = {L0}  T = {d['T']}  pool {d['P']}  E = {d['E']:.3f}  "
              f"est {est_seconds(d['P'], a.threads) / 60:.1f} min", flush=True)
        r = search_multi_fast(L0, factorint(L0), d['T'], threads=a.threads, verbose=True)
        if not r['complete']:
            print("interrupted -- rerun the same command to resume")
            break
        if r['verified'] != r['count']:
            raise RuntimeError(f"unverified hit on L0 = {L0}, T = {d['T']}")
        with open(OUT, 'a') as fh:
            fh.write(json.dumps(dict(L0=str(L0), p0s=r['p0s'], P=r['P'], predicted=r['predicted'], count=r['count'],
                                     verified=r['verified'], seconds=r['seconds'],
                                     values=[str(N) for N, _, _ in r['found']])) + '\n')
        print(f"    found {r['count']} (predicted {r['predicted']:.2f}) in {r['seconds']:.0f}s", flush=True)
        for N, ps, ok in r['found']:
            print(f"    10^{math.log10(N):.2f}  k={len(ps)}  n = {N}")
        runs += 1
    print()
    report()


if __name__ == '__main__':
    main()
