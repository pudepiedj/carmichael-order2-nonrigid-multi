/* cn2kt.c -- memory-lean Wagner k-tree for Howe-construction order-2 Carmichael numbers (driver: cn2_ktree.py --engine c).
 *
 *   cn2kt RUNDIR [threads=32]
 *
 * RUNDIR/job.txt:
 *   K m persist
 *   M_0 ... M_{K-1}                 moduli, pairwise coprime, each < 2^40; M_0..M_{m-1} are the join moduli, the rest filters
 *   2^m lines:  b p_1 ... p_b       leaf i's primes (b <= 32)
 *   2^m lines:  u_0 ... u_{K-1}     leaf i's scale residues (prod over leaves of u == target^-1, so a hit has prod == 1)
 * Output: RUNDIR/hits.txt, one line "H mask_0 ... mask_{2^m-1}" per hit (leaf i's subset); RUNDIR/summary.txt.
 *
 * Tree: node (l, i) at level l joins (l-1, 2i) and (l-1, 2i+1) on a*b == 1 (mod M_{l-1}); the root (level m) keeps
 * the pairs whose filter residues multiply to 1.  Depth-first; every node at levels 1..m-1 is materialised once:
 *   - left leaves are never stored (their residues come from two split tables per modulus, one multiply);
 *   - right leaves are stored as (key, mask), sorted by key;
 *   - a node at level l stores, per entry, its residues modulo M_l .. M_{K-1}, BIT-PACKED (each field exactly
 *     bitlen(M_k) bits, MSB first), sorted by the first field (its next join key), with a bucket index;
 *   - it is built in two passes over its left child (count per bucket, then fill), so it needs no sort scratch.
 * Recovery: for levels >= persist the provenance (index into each child, 5 + 5 bytes) is written to
 * RUNDIR/prov_<l>_<i>.bin in the node's canonical order (sorted by residues, which are distinct); after the run each
 * hit is walked down those files, and the level-1 nodes it needs are rebuilt with their leaf masks (canonical order
 * is reproducible because it depends only on the residues).  The prov files are deleted at the end.
 *
 * build: clang -O3 -mcpu=native -o cn2kt cn2kt.c -lpthread
 */
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>

#define MAXK 16
#define MAXLEAF 256
#define NBITS 20

static int K, m, PERSIST, NT = 32;
static uint64_t MOD[MAXK]; static double INVM[MAXK]; static int W[MAXK];
static char RUNDIR[4096];

typedef struct { int b, h; uint64_t p[32]; uint64_t u[MAXK]; uint64_t *lo[MAXK], *hi[MAXK], *ilo, *ihi; } leaf_t;
static leaf_t LF[MAXLEAF];

static double now(void) { struct timeval t; gettimeofday(&t, 0); return t.tv_sec + 1e-6 * t.tv_usec; }
static double T0;
static double rss_gb(void) { struct rusage r; getrusage(RUSAGE_SELF, &r); return r.ru_maxrss / 1e9; }   /* bytes on macOS */
static void die(const char *s) { fprintf(stderr, "cn2kt: %s\n", s); exit(1); }

/* ---------------- arithmetic (moduli < 2^40) ---------------- */
static inline uint64_t mulmod(uint64_t a, uint64_t b, int k) {
    uint64_t M = MOD[k];
    uint64_t q = (uint64_t)((double)a * (double)b * INVM[k]);
    int64_t r = (int64_t)(a * b - q * M);
    while (r < 0) r += (int64_t)M;
    while (r >= (int64_t)M) r -= (int64_t)M;
    return (uint64_t)r;
}
static uint64_t invmod(uint64_t a, uint64_t M) {
    int64_t t = 0, nt = 1, r = (int64_t)M, nr = (int64_t)(a % M);
    while (nr) { int64_t q = r / nr, x; x = t - q * nt; t = nt; nt = x; x = r - q * nr; r = nr; nr = x; }
    if (r != 1) die("non-invertible residue");
    return (uint64_t)(t < 0 ? t + (int64_t)M : t);
}
static int bitlen(uint64_t x) { int n = 0; while (x) { n++; x >>= 1; } return n; }

/* ---------------- bit-packed fields, MSB first ---------------- */
static inline void putf(uint8_t *r, int o, int w, uint64_t v) {
    while (w > 0) {
        int byte = o >> 3, bo = o & 7, take = 8 - bo; if (take > w) take = w;
        int shift = 8 - bo - take; uint32_t mk = ((1u << take) - 1) << shift;
        r[byte] = (uint8_t)((r[byte] & ~mk) | ((uint32_t)((v >> (w - take)) & ((1u << take) - 1)) << shift));
        o += take; w -= take;
    }
}
static inline uint64_t getf(const uint8_t *r, int o, int w) {
    uint64_t v = 0;
    while (w > 0) {
        int byte = o >> 3, bo = o & 7, take = 8 - bo; if (take > w) take = w;
        int shift = 8 - bo - take;
        v = (v << take) | ((r[byte] >> shift) & ((1u << take) - 1));
        o += take; w -= take;
    }
    return v;
}
static inline void put5(uint8_t *r, uint64_t v) { for (int i = 0; i < 5; i++) r[i] = (uint8_t)(v >> (8 * i)); }
static inline uint64_t get5(const uint8_t *r) { uint64_t v = 0; for (int i = 4; i >= 0; i--) v = v << 8 | r[i]; return v; }

/* ---------------- leaves ---------------- */
static void leaf_tables(leaf_t *L) {
    L->h = L->b / 2;
    int nh = L->b - L->h;
    for (int k = 0; k < K; k++) {
        L->lo[k] = malloc(8 << L->h); L->hi[k] = malloc(8 << nh);
        L->lo[k][0] = L->u[k] % MOD[k]; L->hi[k][0] = 1;
        for (uint64_t s = 1; s < (1ull << L->h); s++) { int t = 63 - __builtin_clzll(s); L->lo[k][s] = mulmod(L->lo[k][s ^ (1ull << t)], L->p[t] % MOD[k], k); }
        for (uint64_t s = 1; s < (1ull << nh); s++) { int t = 63 - __builtin_clzll(s); L->hi[k][s] = mulmod(L->hi[k][s ^ (1ull << t)], L->p[L->h + t] % MOD[k], k); }
    }
    L->ilo = malloc(8 << L->h); L->ihi = malloc(8 << nh);
    L->ilo[0] = invmod(L->u[0] % MOD[0], MOD[0]); L->ihi[0] = 1;
    for (uint64_t s = 1; s < (1ull << L->h); s++) { int t = 63 - __builtin_clzll(s); L->ilo[s] = mulmod(L->ilo[s ^ (1ull << t)], invmod(L->p[t] % MOD[0], MOD[0]), 0); }
    for (uint64_t s = 1; s < (1ull << nh); s++) { int t = 63 - __builtin_clzll(s); L->ihi[s] = mulmod(L->ihi[s ^ (1ull << t)], invmod(L->p[L->h + t] % MOD[0], MOD[0]), 0); }
}
static inline uint64_t leafres(const leaf_t *L, uint64_t mask, int k) {
    return mulmod(L->lo[k][mask & ((1ull << L->h) - 1)], L->hi[k][mask >> L->h], k);
}
static inline uint64_t leafinv0(const leaf_t *L, uint64_t mask) {
    return mulmod(L->ilo[mask & ((1ull << L->h) - 1)], L->ihi[mask >> L->h], 0);
}

/* ---------------- nodes ---------------- */
typedef struct {
    int level, idx;           /* level 0 = a stored right leaf: fields (key mod M_0, mask) */
    int leafid;
    uint64_t n;
    int nf, fk[MAXK + 1], off[MAXK + 1], fw[MAXK + 1];   /* field j: modulus index fk[j] (-1 = mask), bit offset, width */
    int R, Rp;                /* record bytes, and with provenance */
    uint8_t *rec;
    uint64_t *bst; int nbk;   /* bucket starts (2^nbk + 1) over the key field */
} node_t;

static void layout(node_t *X, int level, int prov) {
    X->level = level; X->nf = 0; int o = 0;
    if (level == 0) {
        X->fk[0] = 0; X->fw[0] = W[0]; X->off[0] = 0; o = W[0];
        X->fk[1] = -1; X->fw[1] = 32; X->off[1] = o; o += 32; X->nf = 2;
    } else {
        for (int k = level; k < K; k++) { X->fk[X->nf] = k; X->fw[X->nf] = W[k]; X->off[X->nf] = o; o += W[k]; X->nf++; }
    }
    X->R = (o + 7) / 8; X->Rp = X->R + (prov ? 10 : 0);
    int w0 = X->fw[0]; X->nbk = w0 < NBITS ? w0 : NBITS;
}
static inline uint64_t keyof(const node_t *X, const uint8_t *r) { return getf(r, X->off[0], X->fw[0]); }
static inline uint64_t bucketof(const node_t *X, uint64_t key) { return key >> (X->fw[0] - X->nbk); }

/* field value of modulus k for an entry (node record, or a leaf mask) */
static inline uint64_t nfield(const node_t *X, const uint8_t *r, int k) { return getf(r, X->off[k - X->level], X->fw[k - X->level]); }

static void free_node(node_t *X) { if (!X) return; free(X->rec); free(X->bst); free(X); }

/* sort records of a bucket by their residue bytes (canonical order) */
static int g_cmpR;
static int cmprec(const void *a, const void *b) { return memcmp(a, b, g_cmpR); }
typedef struct { node_t *X; volatile int64_t *next; } sort_arg;
static void *sort_worker(void *v) {
    sort_arg *s = v; node_t *X = s->X; uint64_t nb = 1ull << X->nbk;
    for (;;) {
        int64_t b = __atomic_fetch_add(s->next, 64, __ATOMIC_RELAXED);
        if (b >= (int64_t)nb) break;
        for (int64_t c = b; c < b + 64 && c < (int64_t)nb; c++) {
            uint64_t s0 = X->bst[c], s1 = X->bst[c + 1];
            if (s1 - s0 > 1) qsort(X->rec + s0 * X->Rp, s1 - s0, X->Rp, cmprec);
        }
    }
    return NULL;
}
static void sort_buckets(node_t *X) {
    g_cmpR = X->R; volatile int64_t next = 0; pthread_t th[256]; sort_arg a[256];
    for (int t = 0; t < NT; t++) { a[t] = (sort_arg){X, &next}; pthread_create(&th[t], 0, sort_worker, &a[t]); }
    for (int t = 0; t < NT; t++) pthread_join(th[t], 0);
}

/* ---------------- stored right leaf: (key, mask) sorted by key ---------------- */
typedef struct { node_t *X; const leaf_t *L; int phase; uint32_t *hist; volatile uint64_t *cur; int tid; } la_arg;
static void *leafarr_worker(void *v) {
    la_arg *a = v; node_t *X = a->X; const leaf_t *L = a->L; uint64_t N = 1ull << L->b;
    uint64_t lo = N * a->tid / NT, hi = N * (a->tid + 1) / NT;
    uint32_t *h = a->hist ? a->hist + (uint64_t)a->tid * (1ull << X->nbk) : NULL;
    for (uint64_t s = lo; s < hi; s++) {
        uint64_t key = leafres(L, s, 0), bk = bucketof(X, key);
        if (a->phase == 0) h[bk]++;
        else {
            uint64_t pos = __atomic_fetch_add(&a->cur[bk], 1, __ATOMIC_RELAXED);
            uint8_t *r = X->rec + pos * X->Rp; memset(r, 0, X->Rp);
            putf(r, X->off[0], X->fw[0], key); putf(r, X->off[1], 32, s);
        }
    }
    return NULL;
}
static node_t *build_leafarr(int i) {
    double t = now();
    node_t *X = calloc(1, sizeof *X); layout(X, 0, 0); X->idx = i; X->leafid = i;
    const leaf_t *L = &LF[i]; uint64_t nb = 1ull << X->nbk; X->n = 1ull << L->b;
    uint32_t *hist = calloc((uint64_t)NT * nb, 4);
    pthread_t th[256]; la_arg a[256];
    for (int ph = 0; ph < 2; ph++) {
        volatile uint64_t *cur = NULL;
        if (ph == 1) {
            X->bst = malloc(8 * (nb + 1)); uint64_t s = 0;
            for (uint64_t b = 0; b < nb; b++) { X->bst[b] = s; for (int q = 0; q < NT; q++) s += hist[(uint64_t)q * nb + b]; }
            X->bst[nb] = s;
            X->rec = malloc(X->n * X->Rp + 16);
            cur = malloc(8 * nb); memcpy((void *)cur, X->bst, 8 * nb);
        }
        for (int q = 0; q < NT; q++) { a[q] = (la_arg){X, L, ph, hist, cur, q}; pthread_create(&th[q], 0, leafarr_worker, &a[q]); }
        for (int q = 0; q < NT; q++) pthread_join(th[q], 0);
        if (ph == 1) free((void *)cur);
    }
    free(hist);
    sort_buckets(X);
    fprintf(stderr, "    leaf %d stored: 2^%.2f entries x %d B  (%.0fs, %.0fs total, maxrss %.0f GB)\n", i, (double)L->b, X->Rp,
            now() - t, now() - T0, rss_gb());
    return X;
}

/* ---------------- join ---------------- */
typedef struct {
    int Aleaf; const leaf_t *AL; const node_t *A; uint64_t nA;   /* left: implicit leaf, or a node */
    const node_t *B;                                              /* right: stored leaf (level 0) or node */
    node_t *X;                                                    /* output (NULL at the root) */
    int j;                                                        /* join modulus index (= child level) */
    int phase;                                                    /* 0 count, 1 fill, 2 root */
    uint32_t *hist; volatile uint64_t *cur;
    volatile uint64_t *next;
    pthread_mutex_t *mu; uint64_t *hits; uint64_t *nhits; uint64_t hcap;
    int tid; uint64_t matches;
} join_arg;

static inline uint64_t Bfield(const node_t *B, const uint8_t *rb, int k) {
    if (B->level == 0) return leafres(&LF[B->leafid], getf(rb, B->off[1], 32), k);
    return nfield(B, rb, k);
}
static inline uint64_t Bident(const node_t *B, const uint8_t *rb, uint64_t pos) {
    return B->level == 0 ? getf(rb, B->off[1], 32) : pos;     /* stored leaf: its mask; node: canonical position */
}

#define CH 4096
static void *join_worker(void *v) {
    join_arg *a = v; const node_t *B = a->B; node_t *X = a->X; int j = a->j;
    uint64_t *key = malloc(8 * CH), *pre = malloc(8 * CH);
    uint64_t nbB = 1ull << B->nbk;
    uint32_t *h = (a->phase == 0) ? a->hist + (uint64_t)a->tid * (1ull << X->nbk) : NULL;
    uint64_t matches = 0;
    for (;;) {
        uint64_t s0 = __atomic_fetch_add(a->next, CH, __ATOMIC_RELAXED);
        if (s0 >= a->nA) break;
        uint64_t s1 = s0 + CH < a->nA ? s0 + CH : a->nA, cnt = s1 - s0;
        /* inverse join keys of the left entries */
        if (a->Aleaf) {
            for (uint64_t q = 0; q < cnt; q++) key[q] = leafinv0(a->AL, s0 + q);
        } else {
            uint64_t acc = 1;
            for (uint64_t q = 0; q < cnt; q++) { key[q] = nfield(a->A, a->A->rec + (s0 + q) * a->A->Rp, j); pre[q] = acc; acc = mulmod(acc, key[q], j); }
            uint64_t inv = invmod(acc, MOD[j]);
            for (int64_t q = (int64_t)cnt - 1; q >= 0; q--) { uint64_t kq = key[q]; key[q] = mulmod(inv, pre[q], j); inv = mulmod(inv, kq, j); }
        }
        for (uint64_t q = 0; q < cnt; q++) {
            uint64_t t = key[q], bk = bucketof(B, t);
            if (bk >= nbB) continue;
            uint64_t lo = B->bst[bk], hi = B->bst[bk + 1];
            while (lo < hi) { uint64_t mid = (lo + hi) / 2; if (keyof(B, B->rec + mid * B->Rp) < t) lo = mid + 1; else hi = mid; }
            uint64_t end = B->bst[bk + 1];
            for (uint64_t pb = lo; pb < end; pb++) {
                const uint8_t *rb = B->rec + pb * B->Rp;
                if (keyof(B, rb) != t) break;
                matches++;
                uint64_t ia = s0 + q;
                const uint8_t *ra = a->Aleaf ? NULL : a->A->rec + ia * a->A->Rp;
                #define AF(k) (a->Aleaf ? leafres(a->AL, ia, (k)) : nfield(a->A, ra, (k)))
                if (a->phase == 2) {                       /* root: every filter residue must multiply to 1 */
                    int ok = 1;
                    for (int k = j + 1; k < K && ok; k++) ok = mulmod(AF(k), Bfield(B, rb, k), k) == 1;
                    if (ok) {
                        pthread_mutex_lock(a->mu);
                        if (*a->nhits < a->hcap) { a->hits[2 * *a->nhits] = ia; a->hits[2 * *a->nhits + 1] = Bident(B, rb, pb); (*a->nhits)++; }
                        pthread_mutex_unlock(a->mu);
                    }
                    continue;
                }
                int l = j + 1;
                uint64_t k0 = mulmod(AF(l), Bfield(B, rb, l), l), bx = bucketof(X, k0);
                if (a->phase == 0) { h[bx]++; continue; }
                uint64_t pos = __atomic_fetch_add(&a->cur[bx], 1, __ATOMIC_RELAXED);
                uint8_t *r = X->rec + pos * X->Rp; memset(r, 0, X->Rp);
                putf(r, X->off[0], X->fw[0], k0);
                for (int f = 1; f < X->nf; f++) { int k = X->fk[f]; putf(r, X->off[f], X->fw[f], mulmod(AF(k), Bfield(B, rb, k), k)); }
                if (X->Rp > X->R) { put5(r + X->R, ia); put5(r + X->R + 5, Bident(B, rb, pb)); }
                #undef AF
            }
        }
    }
    a->matches = matches;
    free(key); free(pre);
    return NULL;
}

static uint64_t run_join(join_arg *tmpl, int phase) {
    pthread_t th[256]; join_arg a[256]; volatile uint64_t next = 0; uint64_t tot = 0;
    for (int t = 0; t < NT; t++) { a[t] = *tmpl; a[t].phase = phase; a[t].next = &next; a[t].tid = t; pthread_create(&th[t], 0, join_worker, &a[t]); }
    for (int t = 0; t < NT; t++) { pthread_join(th[t], 0); tot += a[t].matches; }
    return tot;
}

/* materialise node (l, i) from left (implicit leaf or node A) and right B; prov = keep (ia, ib) per entry */
static node_t *join_node(int l, int i, const node_t *A, const node_t *B, int prov) {
    double t = now();
    node_t *X = calloc(1, sizeof *X); layout(X, l, prov); X->idx = i;
    join_arg tm = {0};
    tm.Aleaf = (l == 1); tm.AL = &LF[2 * i]; tm.A = A; tm.nA = (l == 1) ? (1ull << LF[2 * i].b) : A->n;
    tm.B = B; tm.X = X; tm.j = l - 1;
    uint64_t nb = 1ull << X->nbk;
    tm.hist = calloc((uint64_t)NT * nb, 4);
    run_join(&tm, 0);
    X->bst = malloc(8 * (nb + 1)); uint64_t s = 0;
    for (uint64_t b = 0; b < nb; b++) { X->bst[b] = s; for (int q = 0; q < NT; q++) s += tm.hist[(uint64_t)q * nb + b]; }
    X->bst[nb] = s; X->n = s; free(tm.hist); tm.hist = NULL;
    X->rec = malloc(X->n * X->Rp + 16);
    if (!X->rec) die("out of memory for a node");
    volatile uint64_t *cur = malloc(8 * nb); memcpy((void *)cur, X->bst, 8 * nb); tm.cur = cur;
    run_join(&tm, 1);
    free((void *)cur);
    sort_buckets(X);
    fprintf(stderr, "    node (%d,%d): 2^%.2f entries x %d B  (%.0fs, %.0fs total, maxrss %.0f GB)\n", l, i,
            X->n ? __builtin_log2((double)X->n) : 0.0, X->Rp, now() - t, now() - T0, rss_gb());
    return X;
}

/* delete every prov_<l>_<i>.bin of this run directory (levels 2 .. m-1) */
static void clear_prov(void) {
    char fn[4200];
    for (int l = 2; l < m; l++) for (int i = 0; i < (1 << (m - l)); i++) {
        snprintf(fn, sizeof fn, "%s/prov_%d_%d.bin", RUNDIR, l, i); unlink(fn);
    }
}

/* write the provenance of a node to disk in canonical order and drop it from the records */
static void persist_prov(node_t *X) {
    char fn[4200]; snprintf(fn, sizeof fn, "%s/prov_%d_%d.bin", RUNDIR, X->level, X->idx);
    FILE *f = fopen(fn, "wb"); if (!f) die("cannot write a prov file");
    size_t bs = 1 << 20; uint8_t *buf = malloc(10 * bs); size_t nb = 0;
    for (uint64_t e = 0; e < X->n; e++) {
        uint8_t *r = X->rec + e * X->Rp;
        memcpy(buf + 10 * nb, r + X->R, 10); nb++;
        if (nb == bs) { if (fwrite(buf, 10, nb, f) != nb) die("prov write failed (disk full?)"); nb = 0; }
        memmove(X->rec + e * X->R, r, X->R);
    }
    if (nb && fwrite(buf, 10, nb, f) != nb) die("prov write failed (disk full?)");
    fclose(f); free(buf);
    X->Rp = X->R;
    X->rec = realloc(X->rec, X->n * X->R + 16);
}

/* depth-first build of node (l, i), l >= 1; replay = keep provenance in memory (no files) */
static node_t *build(int l, int i, int replay) {
    node_t *A = NULL, *B;
    if (l - 1 >= 1) A = build(l - 1, 2 * i, 0);
    B = (l - 1 == 0) ? build_leafarr(2 * i + 1) : build(l - 1, 2 * i + 1, 0);
    int prov = replay || l >= PERSIST;
    node_t *X = join_node(l, i, A, B, prov);
    free_node(A); free_node(B);
    if (prov && !replay) persist_prov(X);
    return X;
}

static void read_prov(int l, int i, uint64_t e, uint64_t *ia, uint64_t *ib) {
    char fn[4200]; snprintf(fn, sizeof fn, "%s/prov_%d_%d.bin", RUNDIR, l, i);
    int fd = open(fn, O_RDONLY); if (fd < 0) die("cannot read a prov file");
    uint8_t b[10]; if (pread(fd, b, 10, (off_t)(e * 10)) != 10) die("short prov read");
    close(fd); *ia = get5(b); *ib = get5(b + 5);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: cn2kt RUNDIR [threads]\n"); return 1; }
    snprintf(RUNDIR, sizeof RUNDIR, "%s", argv[1]); if (argc > 2) NT = atoi(argv[2]);
    if (NT < 1 || NT > 256) NT = 32;
    char fn[4200]; snprintf(fn, sizeof fn, "%s/job.txt", RUNDIR);
    FILE *f = fopen(fn, "r"); if (!f) die("no job.txt");
    if (fscanf(f, "%d %d %d", &K, &m, &PERSIST) != 3 || K > MAXK || m < 1 || (1 << m) > MAXLEAF) die("bad header");
    for (int k = 0; k < K; k++) {
        if (fscanf(f, "%llu", (unsigned long long *)&MOD[k]) != 1) die("bad modulus");
        if (MOD[k] < 2 || MOD[k] >= (1ull << 40)) die("modulus out of range");
        INVM[k] = 1.0 / (double)MOD[k]; W[k] = bitlen(MOD[k] - 1);
    }
    int nl = 1 << m;
    for (int i = 0; i < nl; i++) {
        if (fscanf(f, "%d", &LF[i].b) != 1 || LF[i].b < 1 || LF[i].b > 32) die("bad leaf");
        for (int q = 0; q < LF[i].b; q++) if (fscanf(f, "%llu", (unsigned long long *)&LF[i].p[q]) != 1) die("bad prime");
    }
    for (int i = 0; i < nl; i++) for (int k = 0; k < K; k++) if (fscanf(f, "%llu", (unsigned long long *)&LF[i].u[k]) != 1) die("bad scale");
    fclose(f);
    T0 = now();
    for (int i = 0; i < nl; i++) leaf_tables(&LF[i]);
    fprintf(stderr, "cn2kt: K=%d m=%d leaves=%d threads=%d persist>=%d\n", K, m, nl, NT, PERSIST);

    /* the run */
    clear_prov();                                  /* stale files of an interrupted run */
    node_t *A = NULL, *B;
    if (m >= 2) A = build(m - 1, 0, 0);
    B = (m == 1) ? build_leafarr(1) : build(m - 1, 1, 0);
    join_arg tm = {0};
    pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER; uint64_t nh = 0, hcap = 1 << 20;
    uint64_t *hits = malloc(16 * hcap);
    tm.Aleaf = (m == 1); tm.AL = &LF[0]; tm.A = A; tm.nA = (m == 1) ? (1ull << LF[0].b) : A->n;
    tm.B = B; tm.j = m - 1; tm.mu = &mu; tm.hits = hits; tm.nhits = &nh; tm.hcap = hcap;
    run_join(&tm, 2);
    free_node(A); free_node(B);
    fprintf(stderr, "cn2kt: root: %llu hits  (%.0fs, maxrss %.0f GB)\n", (unsigned long long)nh, now() - T0, rss_gb());

    /* recovery: masks[h][leaf] */
    uint64_t *masks = calloc(nh * nl + 1, 8);
    if (m == 1) { for (uint64_t h = 0; h < nh; h++) { masks[h * nl] = hits[2 * h]; masks[h * nl + 1] = hits[2 * h + 1]; } }
    else if (nh) {
        /* want[h][node at current level] = entry index; walk the persisted levels */
        uint64_t *want = malloc(8 * nh * nl), *nxt = malloc(8 * nh * nl);
        for (uint64_t h = 0; h < nh; h++) { want[h * nl + 0] = hits[2 * h]; want[h * nl + 1] = hits[2 * h + 1]; }
        int l = m - 1;
        for (; l >= 2 && l >= PERSIST; l--) {
            int nn = 1 << (m - l);
            for (uint64_t h = 0; h < nh; h++) for (int i = 0; i < nn; i++) {
                uint64_t ia, ib; read_prov(l, i, want[h * nl + i], &ia, &ib);
                nxt[h * nl + 2 * i] = ia; nxt[h * nl + 2 * i + 1] = ib;
            }
            uint64_t *t = want; want = nxt; nxt = t;
        }
        if (l != 1) die("persist level must be <= 2 for recovery");
        int nn = 1 << (m - 1);
        for (int i = 0; i < nn; i++) {             /* rebuild each level-1 node with its masks */
            node_t *X = build(1, i, 1);
            for (uint64_t h = 0; h < nh; h++) {
                uint8_t *r = X->rec + want[h * nl + i] * X->Rp;
                masks[h * nl + 2 * i] = get5(r + X->R); masks[h * nl + 2 * i + 1] = get5(r + X->R + 5);
            }
            free_node(X);
        }
        free(want); free(nxt);
    }
    clear_prov();                                  /* hit or null, the provenance (~390 GB at m = 5) is not needed again */
    snprintf(fn, sizeof fn, "%s/hits.txt", RUNDIR);
    f = fopen(fn, "w");
    for (uint64_t h = 0; h < nh; h++) { fprintf(f, "H"); for (int i = 0; i < nl; i++) fprintf(f, " %llu", (unsigned long long)masks[h * nl + i]); fprintf(f, "\n"); }
    fclose(f);
    snprintf(fn, sizeof fn, "%s/summary.txt", RUNDIR);
    f = fopen(fn, "w"); fprintf(f, "hits %llu seconds %.1f maxrss_gb %.1f\n", (unsigned long long)nh, now() - T0, rss_gb()); fclose(f);
    fprintf(stderr, "cn2kt: done, %llu hits, %.0fs\n", (unsigned long long)nh, now() - T0);
    return 0;
}
