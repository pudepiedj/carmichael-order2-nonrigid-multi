# Carmichael numbers of order 2 with two or three exceptional primes

This repository publishes **844 Carmichael numbers of order 2** (OEIS [A175531](https://oeis.org/A175531)) that have *two or three* primes $p$ at which

$$n \equiv p \pmod{p^2-1}$$

instead of $n \equiv 1$. There are **834 with two such primes** and **10 with three**. They were found with Howe's construction, and every one has been verified from scratch.

**These are examples, not a complete list.** A construction finds the numbers on the moduli it is given; it does not sweep everything below a bound. Nothing here says how many such numbers lie below any bound, or that the smallest one found is the smallest that exists.

As far as we know, these are the first substantial lists of such numbers. The first examples, 32 with two such primes and 2 with three, were published with the non-rigid dataset [carmichael-order2-nonrigid](https://github.com/pudepiedj/carmichael-order2-nonrigid) and in the A175531 a-file of non-rigid examples. All 34 are included here.

## The definition

A **Carmichael number of order 2** (Howe 2000) is a squarefree composite $n$ with

$$n \equiv 1 \ \text{ or } \ n \equiv p \pmod{p^2-1} \quad\text{for every prime } p \mid n .$$

| at a prime $p \mid n$ | condition | the map $x \mapsto x^n$ on $\mathbb F_{p^2}$ |
|---|---|---|
| rigid | $n \equiv 1 \pmod{p^2-1}$ | the identity |
| **type $p$** | $n \equiv p \pmod{p^2-1}$ | conjugation, $a+b\sqrt d \mapsto a-b\sqrt d$ |

Below, a **type-$p$ prime** of $n$ means a prime $p \mid n$ with $n \equiv p \pmod{p^2-1}$. A number with at least one is **non-rigid**. Most known non-rigid numbers have exactly one; this repository collects those with two or three. Every number here also has rigid primes: between 21 and 51 of them.

## The smallest examples

| $\log_{10} n$ | primes | type-$p$ primes | prime factors |
|---|---|---|---|
| 66.61 | 24 | 193, 54721 | 31 · 37 · 47 · 79 · 131 · 139 · 181 · **193** · 229 · 233 · 251 · 307 · 463 · 521 · 571 · 599 · 811 · 967 · 2029 · 3191 · 15731 · 19603 · 51679 · **54721** |
| 72.07 | 23 | 1153, 21121 | 53 · 97 · 151 · 311 · 373 · 521 · 811 · 881 · 929 · 991 · **1153** · 1217 · 1301 · 1871 · 2729 · 3191 · 4159 · 4409 · 4523 · 5279 · 5851 · **21121** · 23561 |
| 73.87 | 27 | 2017, 16417 | 43 · 47 · 53 · 61 · 109 · 131 · 139 · 181 · 199 · 239 · 373 · 379 · 433 · 461 · 521 · 571 · 701 · 811 · 919 · 1427 · 1871 · **2017** · 2393 · 2851 · 4003 · **16417** · 62929 |
| 123.29 | 40 | 673, 2017, 2593 | 53 · 61 · 103 · 139 · 181 · 197 · 199 · 239 · 307 · 311 · 379 · 419 · 461 · 463 · 571 · **673** · 701 · 881 · 919 · 967 · 1427 · 1429 · 1451 · 1871 · **2017** · 2393 · 2549 · **2593** · 2729 · 3191 · 4003 · 4049 · 4409 · 4523 · 5851 · 7019 · 11593 · 13339 · 47431 · 244529 |

The first is $n = 4069102870769402987082853488399399844592940278481006187263781896001$. At 22 of its primes $n \equiv 1 \pmod{q^2-1}$; at 193 and 54721, $n \equiv p \pmod{p^2-1}$. Anyone can check this from the factorisation.

## What the dataset holds

| type-$p$ primes per number | numbers | range of $\log_{10} n$ | type-$p$ sets |
|---|---|---|---|
| two | 834 | 66.6 to 148.2 | 86 pairs |
| three | 10 | 123.3 to 160.5 | 4 triples |
| **all** | **844** | | |

The numbers with two type-$p$ primes, by size:

| $\log_{10} n$ | 60s | 70s | 80s | 90s | 100s | 110s | 120s | 130s | 140s |
|---|---|---|---|---|---|---|---|---|---|
| numbers | 1 | 11 | 55 | 151 | 243 | 224 | 113 | 33 | 3 |

The commonest pairs are {193, 6337} (140 numbers), {673, 2017} (76), {673, 2593} (50), {1153, 9601} (47) and {193, 27457} (35). The ten with three type-$p$ primes use {673, 2017, 2593} (6), {673, 2017, 3361} (2), {673, 2017, 7393} (1) and {673, 2593, 8353} (1). These counts reflect where we searched, not how such numbers are distributed.

Every type-$p$ prime here is $\equiv 1 \pmod{24}$. That is necessary whenever $n$ also has a rigid prime, and every number here does.

## Found against predicted

Each search fixes a modulus $L_0$ and a set $T$ of type-$p$ primes. A simple count model predicts its yield as $2^{\lvert P'\rvert}/\varphi(L)$, where $L = \mathrm{lcm}(L_0,\ p^2-1 : p \in T)$ and $P'$ is the pool of primes $q \nmid L_0$ with $q^2 - 1 \mid L_0$, less any prime dividing some $p^2-1$.

| type-$p$ primes | searches | hits found | hits predicted | distinct numbers |
|---|---|---|---|---|
| two | 2,484 | 937 | 900.1 | 834 |
| three | 56 | 10 | 15.1 | 10 |

A hit is one number found by one search. A number found by several searches counts once in the dataset, and its `times_found` says how many searches found it.

The model fits the pairs well. Ten hits against 15.1 predicted for the triples is low, but not alarmingly so: a Poisson count with mean 15.1 is 10 or less about 11% of the time. The triple searches are also not independent, since they share sets $T$ and moduli $L_0$.

## Files

| file | contents |
|---|---|
| `data/cn2_nonrigid_multi.tsv` | the dataset (below) |
| `data/a175531_multi.txt` | the same 844 numbers as an OEIS a-file: the 10 with three type-$p$ primes, then the 834 with two |
| `data/searches/*.jsonl` | the search records: one line per search, with $L_0$, $T$, pool size, prediction, run time and every value found |
| `data/targets/` | the priced sets that were searched: `cn2_twin_pairs.json.gz` (pairs) and `cn2_triple_sets.json` (triples) |
| `data/SHA256SUMS` | SHA-256 hashes of the files above |
| `verify_multi.py` | checks every number and every search record from scratch, standard library only |
| `tools/` | the search programs, as used (below) |

The dataset is tab-separated, with one header line and then one row per number, sorted by $n$:

| column | meaning |
|---|---|
| `log10` | $\log_{10} n$, to 3 decimals |
| `k` | the number of prime factors |
| `type_p` | the type-$p$ primes, comma-separated |
| `L_R` | the rigid part, $\mathrm{lcm}(q^2-1)$ over the rigid primes $q$ |
| `times_found` | how many searches $(L_0, T)$ produced $n$ |
| `n` | the number |
| `primes` | its prime factors, space-separated |

## Checking it

```bash
cd data && shasum -a 256 -c SHA256SUMS && cd ..
python3 verify_multi.py --afile data/a175531_multi.txt --searches data/searches/*.jsonl
```

For every row, `verify_multi.py` checks that:
- the listed factors are distinct primes, and their product is $n$;
- $n \equiv 1$ or $p \pmod{p^2-1}$ at every prime, so $n$ is of order 2;
- the type-$p$ primes are exactly those listed, and there are at least two;
- `L_R` and `log10` are correct;
- the rows are strictly increasing, so no number appears twice.

With `--afile` it checks that every a-file number is in the dataset with the same type-$p$ primes. With `--searches` it reads the search records and checks that:
- every value a search reports is in the dataset, with type-$p$ set $T$ and $n \equiv 1 \pmod{L_0}$;
- `times_found` matches the number of searches reporting each number;
- every number came from some search.

A search that appears twice in the records, because it was rerun, is counted once. The check then prints the found-against-predicted table above. It needs only Python 3 and takes under a second. The last run gave:

```
"rows": 844, "verified": 844, "failed": 0, "order_violations": 0
a-file: 844 numbers, 844 found in the dataset with the same primes p
2 primes p: 2484 searches, 937 hits against 900.1 predicted, 834 distinct numbers
3 primes p: 56 searches, 10 hits against 15.1 predicted, 10 distinct numbers
ALL ROWS VERIFIED
```

## How they were found

Howe's construction fixes a modulus $L_0$ and a type-$p$ set $T$, and looks for

$$n = \Big(\prod_{p \in T} p\Big)\, n_0, \qquad n_0 \text{ squarefree over the pool } P(2,L_0) = \lbrace q \nmid L_0 : q^2-1 \mid L_0\rbrace,$$

with $n \equiv 1 \pmod{L_0}$, $n \equiv 1 \pmod{q^2-1}$ at every $q \mid n_0$, and $n \equiv p \pmod{p^2-1}$ at every $p \in T$. A meet-in-the-middle engine solves each such problem exhaustively, so for a fixed $(L_0, T)$ the search is complete. The dataset is incomplete only because just some pairs $(L_0, T)$ were searched.

The candidates were priced by the count model over the generators of the "box" used for the non-rigid dataset. The cheapest sets per unit of expected yield were then run:

| search | searches | pool sizes $\lvert P'\rvert$ | run time | hits |
|---|---|---|---|---|
| two type-$p$ primes, first searches (September 2026) | 101 | 58–69 | 0.7 h | 33 |
| two type-$p$ primes, harvest (November 2026) | 2,383 | 65–72 | 11.9 h | 904 |
| three type-$p$ primes, first search (September 2026) | 1 | 83 | 8.2 h | 2 |
| three type-$p$ primes, harvest (November 2026) | 56 | 73–84 | 8.9 h | 10 |

The three-prime harvest includes a rerun of the first search, so the two rows share 2 hits.

In the two harvests, every set $T$ was checked against the necessary conditions on sets of type-$p$ primes before it was run:
- every member is $\equiv 1 \pmod{24}$;
- every member has the same 2-adic valuation of $p-1$;
- no odd prime divides both $p+1$ for one member and $p'-1$ for another;
- no member divides another's $p^2-1$.

The programs in `tools/` are the ones used, copied as they were run:

| program | role |
|---|---|
| `cn2x/fast/cn2ss4.c` | the four-list meet-in-the-middle engine in C (compiled on first use with `clang`) |
| `cn2_nonrigid_search_fast.py` | Howe's construction for a set $T$ on one $L_0$, driving `cn2ss4` |
| `cn2_nonrigid_search.py`, `cn2_systematic_fast_SS.py`, `cn2_nonrigid_scout.py` | the original construction, pools and admissibility (the first searches) |
| `cn2_twin_scout.py`, `cn2_multi_scout.py` | price pairs and larger sets on the box generators |
| `cn2_triple_harvest.py` | runs the priced sets, best expected yield per hour first, resumably; checks the conditions above |
| `cn2_nonrigid_collect.py` | merges every search's results, including the one-type-$p$ searches of the earlier dataset, re-verifies them and writes `cn2_nonrigid_all.tsv.gz`; `data/cn2_nonrigid_multi.tsv` is exactly its rows with two or more type-$p$ primes |
| `cn2_multi_afile.py` | writes the a-file from the dataset |

They need `sympy`, `numpy` and `numba` (see `requirements.txt`), and they read and write their records in their own folder. The scouts read the box table `cn2_box_fecundity.csv`, which is not included. Their outputs, the priced sets actually searched, are in `data/targets/`. To rerun a harvest, copy the sets into `tools/`, for example:

```bash
cp data/targets/cn2_triple_sets.json tools/
cd tools
python3 cn2_triple_harvest.py --pmax 84 --hours 8      # resumable: repeat the same command
```

## Licence

MIT. See `LICENSE`.
