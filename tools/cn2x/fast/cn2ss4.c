/*
 * cn2ss4.c -- four-list modular Schroeppel-Shamir for subset products, cache-resident hash join.
 *
 * C successor to the numba kernel of cn2_systematic_fast_SS.py (driver: cn2_nonrigid_search_fast.py).
 * Finds every subset S of a pool P = Q1 u Q2 u Q3 u Q4 with  prod(S) == T  (mod L),  L = A * G * E_1 * ... * E_m
 * pairwise coprime.  For each class t of H = <P mod A>:
 *
 *     list12 = {(i,j) : Q1[i]*Q2[j] == t       (mod A)}   key  Q1[i]*Q2[j]         mod G     (BUILD side)
 *     list34 = {(k,l) : Q3[k]*Q4[l] == T/t     (mod A)}   key  T / (Q3[k]*Q4[l])   mod G     (PROBE side)
 *
 * Equal keys <=> product == T mod A*G; each such candidate is then checked exactly against every modulus
 * (A, G and all E_g) by multiplying the chosen primes out.  Every subset falls in exactly one class, so
 * every solution is found exactly once.
 *
 * What differs from the numba kernel (same mathematics, same answers):
 *   1. UNEQUAL quarters.  Q2 and Q4 are large (up to 2^28 subsets) and bucketed by residue mod A; Q1 and Q3
 *      are small.  A class costs 2^q1 + 2^q3 bucket lookups (cache misses) and each lookup yields a
 *      contiguous run of 2^q2/|H| entries, so the miss count per list entry falls from ~1 to ~1/16.
 *   2. NO SORT.  list12 goes into a bitmap + open-addressing table sized to stay in L2; list34 is never
 *      stored, each entry is generated and probed.  The old kernel built and sorted both lists
 *      (4.5M 64-bit keys each per class at |P| = 83, far outside any cache).
 *   3. Barrett reduction with precomputed constants in place of 64-bit '%'.
 *   4. Per-class lists longer than ~200,000 entries (pools above ~90) are hash-partitioned before the join, so
 *      the join tables stay in cache at any list length (CN2SS4_PARTMIN overrides the threshold).
 *   5. Resumable: completed class-chunks are logged; rerun the same command to continue.
 *
 * Usage:  cn2ss4 RUNDIR [threads=32] [progress_sec=30]
 *   RUNDIR/job.txt   n q1 q2 q3 q4 / primes / A TA / G TG / ne / E_g TE_g ...        (written by the driver)
 *   RUNDIR/out.txt   "H i j k l"  -- subset masks in Q1..Q4 of each solution (may repeat after a resume)
 *   RUNDIR/done.log  "C <chunk>"  per finished chunk (~2^31 list entries each);   RUNDIR/STOP or SIGINT/SIGTERM stops cleanly
 *   RUNDIR/summary.txt  written only when every chunk is done
 *
 * Build:  clang -O3 -mcpu=native -o cn2ss4 cn2ss4.c -lpthread -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <stdatomic.h>

typedef uint32_t u32;
typedef uint64_t u64;
typedef unsigned __int128 u128;

#define MAXN 128
#define MAXE 8
#define BLK 32                       /* lookups issued together so their cache misses overlap */

static void die(const char *m) { fprintf(stderr, "cn2ss4: %s\n", m); exit(1); }
static void *xmalloc(size_t n) { void *p = malloc(n ? n : 1); if (!p) die("out of memory"); return p; }
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

/* ---- Barrett: x < 2^64, 1 <= M < 2^32, mu = floor((2^64-1)/M);  q >= floor(x/M) - 1 so one correction ---- */
static inline u64 bred(u64 x, u64 M, u64 mu) { u64 r = x - (u64)(((u128)x * mu) >> 64) * M; return r >= M ? r - M : r; }
static inline u64 mulmod64(u64 a, u64 b, u64 M) { return (u64)((u128)a * b % M); }
static u64 modinv(u64 a, u64 m) {
    if (m == 1) return 0;
    __int128 r0 = m, r1 = a % m, x0 = 0, x1 = 1;
    while (r1) { __int128 q = r0 / r1, t = r0 - q * r1; r0 = r1; r1 = t; t = x0 - q * x1; x0 = x1; x1 = t; }
    if (r0 != 1) die("modinv: not coprime");
    return (u64)(x0 < 0 ? x0 + m : x0);
}

/* ---- job ---- */
static int n, qs[4], qoff[4], ne;
static u64 P[MAXN], A, TA, G, TG, E[MAXE], TE[MAXE], muA, muG;
static u32 N1, N2, N3, N4;

/* ---- tables ---- */
static u32 *ia1, *g1;                /* Q1: (prod)^-1 mod A,  prod mod G */
static u32 *ia3, *ig3;               /* Q3: (prod)^-1 mod A,  TG * (prod)^-1 mod G */
static u32 *off2, *g2o, *ord2;       /* Q2 bucketed by prod mod A: offsets, prod mod G, subset mask */
static u32 *off4, *ig4o, *ord4;      /* Q4 bucketed by prod mod A: offsets, (prod)^-1 mod G, subset mask */
static u32 *H; static u64 h;         /* classes */

static void subset_table(u32 *out, const u64 *pr, int q, u64 M, int inverse) {
    u64 mu = M > 1 ? (~(u64)0) / M : ~(u64)0;
    out[0] = (u32)(1 % M);
    for (int b = 0; b < q; b++) {
        u64 pb = inverse ? modinv(pr[b] % M, M) : pr[b] % M;
        size_t half = (size_t)1 << b;
        for (size_t m = 0; m < half; m++) out[half + m] = (u32)bred(out[m] * pb, M, mu);
    }
}

typedef struct { const u64 *pr; int q, ginv; u32 **off, **go, **ord; } bucket_job;
static void *bucket_thread(void *arg) {                 /* counting sort of the 2^q subsets by residue mod A */
    bucket_job *j = arg;
    size_t N = (size_t)1 << j->q;
    u32 *rA = xmalloc(N * 4), *g = xmalloc(N * 4);
    subset_table(rA, j->pr, j->q, A, 0);
    subset_table(g, j->pr, j->q, G, j->ginv);
    u32 *off = calloc(A + 2, 4); if (!off) die("out of memory");
    for (size_t i = 0; i < N; i++) off[rA[i] + 1]++;
    for (u64 r = 0; r < A; r++) off[r + 1] += off[r];
    u32 *pos = xmalloc(A * 4); memcpy(pos, off, A * 4);
    u32 *go = xmalloc(N * 4), *ord = xmalloc(N * 4);
    for (size_t i = 0; i < N; i++) { u32 p = pos[rA[i]]++; go[p] = g[i]; ord[p] = (u32)i; }
    free(pos); free(rA); free(g);
    *j->off = off; *j->go = go; *j->ord = ord;
    return NULL;
}

static void build_classes(void) {                       /* BFS of <P mod A> */
    unsigned char *seen = calloc(A, 1); if (!seen) die("out of memory");
    u32 gens[MAXN]; int ng = 0;
    for (int i = 0; i < n; i++) { u32 g = (u32)(P[i] % A); int dup = 0; for (int k = 0; k < ng; k++) dup |= gens[k] == g; if (!dup) gens[ng++] = g; }
    size_t capq = 1 << 16; H = xmalloc(capq * 4);
    u32 s = (u32)(1 % A); seen[s] = 1; H[0] = s; size_t head = 0, tail = 1;
    while (head < tail) {
        u64 x = H[head++];
        for (int k = 0; k < ng; k++) {
            u32 y = (u32)bred(x * gens[k], A, muA);
            if (!seen[y]) { seen[y] = 1; if (tail == capq) { capq *= 2; H = realloc(H, capq * 4); if (!H) die("out of memory"); } H[tail++] = y; }
        }
    }
    free(seen); h = tail;
}

/* ---- run state ---- */
static const char *rundir;
static u64 nchunks, chunk_classes;
static unsigned char *chunk_done;
static atomic_ullong next_chunk, chunks_finished, tot_cand, tot_hits, tot_e12, tot_e34, tot_maxlist;
static atomic_int stop_flag;
static pthread_mutex_t io_mu = PTHREAD_MUTEX_INITIALIZER;
static FILE *f_out, *f_done;
static u64 mean12;

static void on_sig(int s) { (void)s; atomic_store(&stop_flag, 1); }

static int verify(u32 i, u32 j, u32 k, u32 l) {         /* exact: every modulus, primes multiplied out */
    u32 mk[4] = { i, j, k, l };
    u64 mods[MAXE + 2], tg[MAXE + 2]; int nm = 0;
    mods[nm] = A; tg[nm++] = TA; mods[nm] = G; tg[nm++] = TG;
    for (int g = 0; g < ne; g++) { mods[nm] = E[g]; tg[nm++] = TE[g]; }
    for (int m = nm - 1; m >= 0; m--) {                 /* E first: they are the untested ones */
        u64 M = mods[m], v = 1 % M;
        for (int q = 0; q < 4; q++) for (int b = 0; b < qs[q]; b++) if (mk[q] >> b & 1) v = mulmod64(v, P[qoff[q] + b] % M, M);
        if (v != tg[m]) { if (m < 2) die("internal error: join candidate fails A or G"); return 0; }
    }
    return 1;
}

typedef struct { u64 *bm; u32 *T, *K, *I, *U; u64 bm_words, slots, cap; } wspace;
static void ws_alloc(wspace *w, u64 want) {
    u64 slots = 64; while (slots < 2 * want) slots <<= 1;
    w->slots = slots; w->cap = slots - slots / 4; w->bm_words = slots / 8;      /* bitmap = 8 bits per slot */
    w->bm = xmalloc(w->bm_words * 8); w->T = xmalloc(slots * 4);
    w->K = xmalloc(w->cap * 4); w->I = xmalloc(w->cap * 4); w->U = xmalloc(w->cap * 4);
}
static void ws_free(wspace *w) { free(w->bm); free(w->T); free(w->K); free(w->I); free(w->U); }

#define HASH(k) ((u32)(k) * 0x9E3779B1u)

/* one class; returns 0 if list12 overflowed the workspace (caller grows and retries) */
static int do_class(wspace *w, u64 t, u64 *cand, u64 *e12, u64 *e34, u64 *hits) {
    const u64 a = A, ma = muA, g = G, mg = muG;
    u64 tinv = bred(TA * modinv(t, a), a, ma);
    memset(w->bm, 0, w->bm_words * 8);
    memset(w->T, 0, w->slots * 4);
    const int tsh = 32 - __builtin_ctzll(w->slots), bsh = 32 - __builtin_ctzll(w->bm_words * 64);
    const u32 smask = (u32)(w->slots - 1);
    u64 *bm = w->bm; u32 *T = w->T, *K = w->K, *I = w->I, *U = w->U;
    u64 m = 0; const u64 cap = w->cap;
    u32 rr[BLK], lo[BLK], hi[BLK];

    for (u32 i0 = 0; i0 < N1; i0 += BLK) {                              /* BUILD list12 */
        u32 nb = N1 - i0 < BLK ? N1 - i0 : BLK;
        for (u32 x = 0; x < nb; x++) { rr[x] = (u32)bred(t * ia1[i0 + x], a, ma); __builtin_prefetch(&off2[rr[x]]); }
        for (u32 x = 0; x < nb; x++) { lo[x] = off2[rr[x]]; hi[x] = off2[rr[x] + 1]; __builtin_prefetch(&g2o[lo[x]]); }
        for (u32 x = 0; x < nb; x++) {
            u64 gi = g1[i0 + x];
            if (m + (hi[x] - lo[x]) > cap) return 0;
            for (u32 u = lo[x]; u < hi[x]; u++) {
                u32 key = (u32)bred(gi * g2o[u], g, mg), hh = HASH(key);
                bm[(hh >> bsh) >> 6] |= (u64)1 << ((hh >> bsh) & 63);
                u32 s = hh >> tsh; while (T[s]) s = (s + 1) & smask;
                T[s] = (u32)m + 1; K[m] = key; I[m] = i0 + x; U[m] = u; m++;
            }
        }
    }
    *e12 += m;
    if (m == 0) return 1;

    u64 n34 = 0;
    for (u32 k0 = 0; k0 < N3; k0 += BLK) {                              /* PROBE with list34 */
        u32 nb = N3 - k0 < BLK ? N3 - k0 : BLK;
        for (u32 x = 0; x < nb; x++) { rr[x] = (u32)bred(tinv * ia3[k0 + x], a, ma); __builtin_prefetch(&off4[rr[x]]); }
        for (u32 x = 0; x < nb; x++) { lo[x] = off4[rr[x]]; hi[x] = off4[rr[x] + 1]; __builtin_prefetch(&ig4o[lo[x]]); }
        for (u32 x = 0; x < nb; x++) {
            u64 gk = ig3[k0 + x];
            n34 += hi[x] - lo[x];
            for (u32 u = lo[x]; u < hi[x]; u++) {
                u32 key = (u32)bred(gk * ig4o[u], g, mg), hh = HASH(key);
                if (!(bm[(hh >> bsh) >> 6] >> ((hh >> bsh) & 63) & 1)) continue;
                for (u32 s = hh >> tsh; T[s]; s = (s + 1) & smask) {
                    u32 e = T[s] - 1;
                    if (K[e] != key) continue;
                    u32 i = I[e], j = ord2[U[e]], kk = k0 + x, l = ord4[u];
                    if (!(i | j | kk | l)) continue;                    /* empty product */
                    ++*cand;
                    if (verify(i, j, kk, l)) {
                        ++*hits;
                        pthread_mutex_lock(&io_mu);
                        fprintf(f_out, "H %u %u %u %u\n", i, j, kk, l); fflush(f_out);
                        pthread_mutex_unlock(&io_mu);
                    }
                }
            }
        }
    }
    *e34 += n34;
    return 1;
}


/* ---- partitioned join, for per-class lists too long for L2 ----
 * Both lists are generated straight into 2^kbits partitions by the top bits of the key hash (sequential write
 * streams), then each partition is joined through a table of ~2^15 slots.  Cost per entry is flat in the list length. */
typedef struct { u32 key, a, b; } rec;
typedef struct { rec *L12, *L34; u32 *n12, *n34, *T; u64 cap12, cap34, slots; int kbits; } pspace;
static void ps_alloc(pspace *w, u64 m12, u64 m34, int grow) {
    int kb = 1; while (kb < 12 && (m12 >> kb) > 8192) kb++;
    w->kbits = kb;
    w->cap12 = ((m12 >> kb) + 256) * 3 / 2 << grow; w->cap34 = ((m34 >> kb) + 256) * 3 / 2 << grow;
    u64 slots = 64; while (slots < 2 * w->cap12) slots <<= 1;
    w->slots = slots;
    w->L12 = xmalloc((w->cap12 << kb) * sizeof(rec)); w->L34 = xmalloc((w->cap34 << kb) * sizeof(rec));
    w->n12 = xmalloc(4u << kb); w->n34 = xmalloc(4u << kb); w->T = xmalloc(slots * 4);
}
static void ps_free(pspace *w) { free(w->L12); free(w->L34); free(w->n12); free(w->n34); free(w->T); }

static int do_class_part(pspace *w, u64 t, u64 *cand, u64 *e12, u64 *e34, u64 *hits) {
    const u64 a = A, ma = muA, g = G, mg = muG;
    u64 tinv = bred(TA * modinv(t, a), a, ma);
    const int kb = w->kbits, psh = 32 - kb;
    const u64 c12 = w->cap12, c34 = w->cap34;
    u32 *n12 = w->n12, *n34 = w->n34, *T = w->T;
    memset(n12, 0, 4u << kb); memset(n34, 0, 4u << kb);
    u32 rr[BLK], lo[BLK], hi[BLK];
    u64 m12 = 0, m34 = 0;
    for (u32 i0 = 0; i0 < N1; i0 += BLK) {
        u32 nb = N1 - i0 < BLK ? N1 - i0 : BLK;
        for (u32 x = 0; x < nb; x++) { rr[x] = (u32)bred(t * ia1[i0 + x], a, ma); __builtin_prefetch(&off2[rr[x]]); }
        for (u32 x = 0; x < nb; x++) { lo[x] = off2[rr[x]]; hi[x] = off2[rr[x] + 1]; __builtin_prefetch(&g2o[lo[x]]); }
        for (u32 x = 0; x < nb; x++) {
            u64 gi = g1[i0 + x];
            for (u32 u = lo[x]; u < hi[x]; u++) {
                u32 key = (u32)bred(gi * g2o[u], g, mg), p = HASH(key) >> psh;
                if (n12[p] == c12) return 0;
                rec *r = &w->L12[p * c12 + n12[p]++]; r->key = key; r->a = i0 + x; r->b = u;
            }
            m12 += hi[x] - lo[x];
        }
    }
    for (u32 k0 = 0; k0 < N3; k0 += BLK) {
        u32 nb = N3 - k0 < BLK ? N3 - k0 : BLK;
        for (u32 x = 0; x < nb; x++) { rr[x] = (u32)bred(tinv * ia3[k0 + x], a, ma); __builtin_prefetch(&off4[rr[x]]); }
        for (u32 x = 0; x < nb; x++) { lo[x] = off4[rr[x]]; hi[x] = off4[rr[x] + 1]; __builtin_prefetch(&ig4o[lo[x]]); }
        for (u32 x = 0; x < nb; x++) {
            u64 gk = ig3[k0 + x];
            for (u32 u = lo[x]; u < hi[x]; u++) {
                u32 key = (u32)bred(gk * ig4o[u], g, mg), p = HASH(key) >> psh;
                if (n34[p] == c34) return 0;
                rec *r = &w->L34[p * c34 + n34[p]++]; r->key = key; r->a = k0 + x; r->b = u;
            }
            m34 += hi[x] - lo[x];
        }
    }
    const u32 smask = (u32)(w->slots - 1);
    for (u32 p = 0; p < (1u << kb); p++) {
        u32 na = n12[p], nb = n34[p];
        if (!na || !nb) continue;
        const rec *RA = &w->L12[p * c12], *RB = &w->L34[p * c34];
        memset(T, 0, w->slots * 4);
        for (u32 e = 0; e < na; e++) { u32 s = (HASH(RA[e].key) >> 4) & smask; while (T[s]) s = (s + 1) & smask; T[s] = e + 1; }
        for (u32 y = 0; y < nb; y++) {
            u32 key = RB[y].key;
            for (u32 s = (HASH(key) >> 4) & smask; T[s]; s = (s + 1) & smask) {
                const rec *ra = &RA[T[s] - 1];
                if (ra->key != key) continue;
                u32 i = ra->a, j = ord2[ra->b], kk = RB[y].a, l = ord4[RB[y].b];
                if (!(i | j | kk | l)) continue;
                ++*cand;
                if (verify(i, j, kk, l)) {
                    ++*hits;
                    pthread_mutex_lock(&io_mu);
                    fprintf(f_out, "H %u %u %u %u\n", i, j, kk, l); fflush(f_out);
                    pthread_mutex_unlock(&io_mu);
                }
            }
        }
    }
    *e12 += m12; *e34 += m34;
    return 1;
}

static int use_part; static u64 mean34;
static void *worker(void *arg) {
    (void)arg;
    wspace w; pspace pw; int grow = 0;
    if (use_part) ps_alloc(&pw, mean12, mean34, 0); else ws_alloc(&w, mean12 + mean12 / 4 + 64);
    for (;;) {
        u64 c = atomic_fetch_add(&next_chunk, 1);
        if (c >= nchunks || atomic_load(&stop_flag)) break;
        if (chunk_done[c]) continue;
        u64 lo = c * chunk_classes, hi = lo + chunk_classes; if (hi > h) hi = h;
        u64 cand = 0, e12 = 0, e34 = 0, hits = 0, maxl = 0; int aborted = 0;
        for (u64 ci = lo; ci < hi; ci++) {
            if (atomic_load(&stop_flag)) { aborted = 1; break; }
            u64 before = e12;
            if (use_part) {
                u64 c0 = cand, h0 = hits;                               /* a partition overflowed: double and redo */
                while (!do_class_part(&pw, H[ci], &cand, &e12, &e34, &hits)) { cand = c0; hits = h0; ps_free(&pw); ps_alloc(&pw, mean12, mean34, ++grow); }
            } else while (!do_class(&w, H[ci], &cand, &e12, &e34, &hits)) {    /* list12 too long: double and redo */
                u64 want = w.slots; ws_free(&w); ws_alloc(&w, want);
            }
            if (e12 - before > maxl) maxl = e12 - before;
        }
        if (aborted) break;                                             /* chunk not logged: redone on resume */
        pthread_mutex_lock(&io_mu);
        fprintf(f_done, "C %llu\n", (unsigned long long)c); fflush(f_done);
        pthread_mutex_unlock(&io_mu);
        atomic_fetch_add(&tot_cand, cand); atomic_fetch_add(&tot_hits, hits);
        atomic_fetch_add(&tot_e12, e12); atomic_fetch_add(&tot_e34, e34);
        u64 pm = atomic_load(&tot_maxlist); while (maxl > pm && !atomic_compare_exchange_weak(&tot_maxlist, &pm, maxl)) ;
        atomic_fetch_add(&chunks_finished, 1);
    }
    if (use_part) ps_free(&pw); else ws_free(&w);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 2) die("usage: cn2ss4 RUNDIR [threads=32] [progress_sec=30]");
    rundir = argv[1];
    int threads = argc > 2 ? atoi(argv[2]) : 32;
    double prog = argc > 3 ? atof(argv[3]) : 30.0;
    char path[4096];
    snprintf(path, sizeof path, "%s/job.txt", rundir);
    FILE *f = fopen(path, "r"); if (!f) die("cannot open RUNDIR/job.txt");
    if (fscanf(f, "%d %d %d %d %d", &n, &qs[0], &qs[1], &qs[2], &qs[3]) != 5) die("bad job header");
    if (n < 4 || n > MAXN || qs[0] + qs[1] + qs[2] + qs[3] != n) die("bad n / quarters");
    for (int q = 0; q < 4; q++) if (qs[q] < 1 || qs[q] > 31) die("quarter size out of range 1..31");
    for (int i = 0; i < n; i++) if (fscanf(f, "%llu", (unsigned long long *)&P[i]) != 1) die("bad primes");
    unsigned long long a_, ta_, g_, tg_;
    if (fscanf(f, "%llu %llu %llu %llu %d", &a_, &ta_, &g_, &tg_, &ne) != 5) die("bad A/G lines");
    A = a_; TA = ta_; G = g_; TG = tg_;
    if (A < 1 || A >= ((u64)1 << 31) || G < 1 || G >= ((u64)1 << 32) || ne < 0 || ne > MAXE) die("A, G or ne out of range");
    for (int g = 0; g < ne; g++) { unsigned long long e_, te_; if (fscanf(f, "%llu %llu", &e_, &te_) != 2) die("bad E line"); E[g] = e_; TE[g] = te_; if (E[g] >> 63) die("E too large"); }
    fclose(f);
    qoff[0] = 0; for (int q = 1; q < 4; q++) qoff[q] = qoff[q - 1] + qs[q - 1];
    N1 = 1u << qs[0]; N2 = 1u << qs[1]; N3 = 1u << qs[2]; N4 = 1u << qs[3];
    muA = A > 1 ? (~(u64)0) / A : ~(u64)0; muG = G > 1 ? (~(u64)0) / G : ~(u64)0;

    double t0 = now();
    ia1 = xmalloc((size_t)N1 * 4); g1 = xmalloc((size_t)N1 * 4); ia3 = xmalloc((size_t)N3 * 4); ig3 = xmalloc((size_t)N3 * 4);
    subset_table(ia1, P + qoff[0], qs[0], A, 1); subset_table(g1, P + qoff[0], qs[0], G, 0);
    subset_table(ia3, P + qoff[2], qs[2], A, 1); subset_table(ig3, P + qoff[2], qs[2], G, 1);
    for (u32 k = 0; k < N3; k++) ig3[k] = (u32)bred((u64)ig3[k] * TG, G, muG);
    bucket_job b2 = { P + qoff[1], qs[1], 0, &off2, &g2o, &ord2 }, b4 = { P + qoff[3], qs[3], 1, &off4, &ig4o, &ord4 };
    pthread_t tb; pthread_create(&tb, NULL, bucket_thread, &b2); bucket_thread(&b4); pthread_join(tb, NULL);
    build_classes();
    mean12 = (u64)ceil(ldexp(1.0, qs[0] + qs[1]) / (double)h);
    mean34 = (u64)ceil(ldexp(1.0, qs[2] + qs[3]) / (double)h);
    { const char *e = getenv("CN2SS4_PARTMIN"); use_part = mean12 > (u64)(e ? atof(e) : 200000.0); }
    fprintf(stderr, "cn2ss4: |P| = %d  quarters %d %d %d %d  A = %llu  |H| = %llu (2^%.2f)  G = %llu  extra moduli %d\n"
            "        list12 ~%llu per class, list34 ~%.0f (%s join); tables built in %.1fs\n",
            n, qs[0], qs[1], qs[2], qs[3], a_, (unsigned long long)h, log2((double)h), g_, ne,
            (unsigned long long)mean12, ldexp(1.0, qs[2] + qs[3]) / (double)h, use_part ? "partitioned" : "direct", now() - t0);

    /* a chunk is ~2^31 list entries (10-30 s of one thread) and never more than 1/1024 of the classes; independent of threads */
    chunk_classes = (u64)(ldexp(1.0, 31) / (ldexp(1.0, qs[0] + qs[1]) / (double)h + ldexp(1.0, qs[2] + qs[3]) / (double)h));
    if (chunk_classes > h / 1024) chunk_classes = h / 1024;
    if (chunk_classes < 1) chunk_classes = 1;
    nchunks = (h + chunk_classes - 1) / chunk_classes;
    chunk_done = calloc(nchunks, 1);
    u64 already = 0;
    snprintf(path, sizeof path, "%s/done.log", rundir);
    if ((f = fopen(path, "r"))) {
        unsigned long long cc, hh, c;
        if (fscanf(f, "J %llu %llu", &hh, &cc) == 2) {
            if (hh != h || cc != chunk_classes) die("done.log was written with a different |H| or chunk size; use a new RUNDIR");
            while (fscanf(f, " C %llu", &c) == 1) if (c < nchunks && !chunk_done[c]) { chunk_done[c] = 1; already++; }
        }
        fclose(f);
    }
    f_done = fopen(path, "a"); if (!f_done) die("cannot open done.log");
    if (already == 0 && ftell(f_done) == 0) { fprintf(f_done, "J %llu %llu\n", (unsigned long long)h, (unsigned long long)chunk_classes); fflush(f_done); }
    snprintf(path, sizeof path, "%s/out.txt", rundir);
    f_out = fopen(path, "a"); if (!f_out) die("cannot open out.txt");
    if (already) fprintf(stderr, "cn2ss4: resuming, %llu of %llu chunks already done\n", (unsigned long long)already, (unsigned long long)nchunks);

    signal(SIGINT, on_sig); signal(SIGTERM, on_sig);
    pthread_t *th = xmalloc(sizeof(pthread_t) * threads);
    double t1 = now(), tlast = t1;
    for (int i = 0; i < threads; i++) pthread_create(&th[i], NULL, worker, NULL);
    char stop_path[4096]; snprintf(stop_path, sizeof stop_path, "%s/STOP", rundir);
    u64 todo = nchunks - already;
    for (;;) {
        usleep(200000);
        u64 fin = atomic_load(&chunks_finished);
        if (access(stop_path, F_OK) == 0) atomic_store(&stop_flag, 1);
        if (fin >= todo || atomic_load(&stop_flag)) break;
        if (now() - tlast >= prog) {
            tlast = now(); double el = tlast - t1, fr = todo ? (double)fin / todo : 1;
            fprintf(stderr, "  %6.2f%%  %llu hits  %llu candidates  %.0fs  ETA %.0fs\n", 100.0 * (already + fin) / nchunks,
                    (unsigned long long)atomic_load(&tot_hits), (unsigned long long)atomic_load(&tot_cand), el, fr > 0 ? el / fr - el : 0.0);
        }
    }
    for (int i = 0; i < threads; i++) pthread_join(th[i], NULL);
    double el = now() - t1;
    u64 fin = atomic_load(&chunks_finished);
    int complete = already + fin == nchunks;
    fprintf(stderr, "cn2ss4: %s  %llu/%llu chunks  this session: %llu hits, %llu candidates, %.3g list entries, max list12 %llu, %.1fs (%.2f ns/entry/thread)\n",
            complete ? "COMPLETE" : "STOPPED (rerun to resume)", (unsigned long long)(already + fin), (unsigned long long)nchunks,
            (unsigned long long)atomic_load(&tot_hits), (unsigned long long)atomic_load(&tot_cand),
            (double)(atomic_load(&tot_e12) + atomic_load(&tot_e34)), (unsigned long long)atomic_load(&tot_maxlist), el,
            1e9 * el * threads / (double)(atomic_load(&tot_e12) + atomic_load(&tot_e34) + 1));
    if (complete) {
        snprintf(path, sizeof path, "%s/summary.txt", rundir);
        f = fopen(path, "a");
        fprintf(f, "complete 1\nH %llu\nA %llu\nG %llu\nsession_seconds %.1f\nsession_candidates %llu\nsession_entries %.6g\nthreads %d\n",
                (unsigned long long)h, a_, g_, el + (t1 - t0), (unsigned long long)atomic_load(&tot_cand),
                (double)(atomic_load(&tot_e12) + atomic_load(&tot_e34)), threads);
        fclose(f);
    }
    fclose(f_out); fclose(f_done);
    return complete ? 0 : 2;
}
