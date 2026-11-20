/* ktjoin.c -- threaded equi-join of two uint32 key arrays, for cn2_ktree.py (Wagner k-tree).
 *
 *   int64_t kt_join(const uint32_t *ka, int64_t nA, const uint32_t *kb, int64_t nB,
 *                   uint32_t *outA, uint32_t *outB, int64_t cap, int threads)
 * writes every pair (i, j) with ka[i] == kb[j] into outA/outB and returns the number of pairs;
 * if that exceeds cap nothing is guaranteed in out*, and the caller retries with a larger cap.
 *
 *   void kt_mulgather(const uint32_t *x, const uint32_t *y, const uint32_t *ia, const uint32_t *ib,
 *                     int64_t n, uint32_t M, uint32_t *out, int threads)
 * out[k] = x[ia[k]] * y[ib[k]] mod M.
 *
 * Method: both sides are radix-partitioned on the low key bits into buckets of a few thousand entries
 * (parallel counting sort), then each bucket is joined with a small hash table that stays in cache.
 * Two passes over the buckets: count, prefix sum, fill.
 *
 * build: clang -O3 -mcpu=native -shared -fPIC -o libktjoin.dylib ktjoin.c -lpthread
 */
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const uint32_t *keys; int64_t n; int pbits; int nt;
    int64_t *hist;          /* nt * nb */
    uint32_t *okey, *oidx;
    int tid; int phase;
} part_t;

static void *part_worker(void *v) {
    part_t *p = v;
    int64_t nb = (int64_t)1 << p->pbits, mask = nb - 1;
    int64_t lo = p->n * p->tid / p->nt, hi = p->n * (p->tid + 1) / p->nt;
    int64_t *h = p->hist + (int64_t)p->tid * nb;
    if (p->phase == 0) {
        memset(h, 0, nb * sizeof(int64_t));
        for (int64_t i = lo; i < hi; i++) h[p->keys[i] & mask]++;
    } else {
        for (int64_t i = lo; i < hi; i++) {
            uint32_t k = p->keys[i];
            int64_t d = h[k & mask]++;
            p->okey[d] = k; p->oidx[d] = (uint32_t)i;
        }
    }
    return NULL;
}

/* partition keys into nb buckets; off has nb+1 entries */
static void partition(const uint32_t *keys, int64_t n, int pbits, int nt,
                      uint32_t *okey, uint32_t *oidx, int64_t *off) {
    int64_t nb = (int64_t)1 << pbits;
    int64_t *hist = malloc(sizeof(int64_t) * nb * nt);
    pthread_t th[256]; part_t a[256];
    for (int ph = 0; ph < 2; ph++) {
        for (int t = 0; t < nt; t++) {
            a[t] = (part_t){keys, n, pbits, nt, hist, okey, oidx, t, ph};
            pthread_create(&th[t], NULL, part_worker, &a[t]);
        }
        for (int t = 0; t < nt; t++) pthread_join(th[t], NULL);
        if (ph == 0) {          /* column-major prefix: bucket b, thread t */
            int64_t s = 0;
            for (int64_t b = 0; b < nb; b++) {
                off[b] = s;
                for (int t = 0; t < nt; t++) { int64_t c = hist[(int64_t)t * nb + b]; hist[(int64_t)t * nb + b] = s; s += c; }
            }
            off[nb] = s;
        }
    }
    free(hist);
}

typedef struct {
    const uint32_t *akey, *aidx, *bkey, *bidx; const int64_t *offA, *offB;
    int64_t nb; int pbits; int phase;
    int64_t *bcnt;          /* per bucket: pairs (phase 0), then output offset (phase 1) */
    uint32_t *outA, *outB;
    volatile int64_t *next;
} join_t;

static inline uint32_t hsh(uint32_t k, int pbits) { return (k >> pbits) * 2654435761u; }

static void *join_worker(void *v) {
    join_t *J = v;
    int64_t cap = 0; int32_t *tab = NULL, *chain = NULL; int64_t ccap = 0;
    for (;;) {
        int64_t b = __atomic_fetch_add(J->next, 1, __ATOMIC_RELAXED);
        if (b >= J->nb) break;
        int64_t a0 = J->offA[b], a1 = J->offA[b + 1], b0 = J->offB[b], b1 = J->offB[b + 1];
        int64_t nB = b1 - b0;
        if (nB == 0 || a1 == a0) { if (J->phase == 0) J->bcnt[b] = 0; continue; }
        int64_t sz = 16; while (sz < 2 * nB) sz <<= 1;
        if (sz > cap) { free(tab); cap = sz; tab = malloc(sizeof(int32_t) * cap); }
        if (nB > ccap) { free(chain); ccap = nB; chain = malloc(sizeof(int32_t) * ccap); }
        memset(tab, 0xff, sizeof(int32_t) * sz);
        uint32_t m = (uint32_t)(sz - 1);
        /* table holds the head of a chain of equal keys */
        for (int64_t j = 0; j < nB; j++) {
            uint32_t k = J->bkey[b0 + j], h = hsh(k, J->pbits) & m;
            for (;;) {
                int32_t e = tab[h];
                if (e < 0) { tab[h] = (int32_t)j; chain[j] = -1; break; }
                if (J->bkey[b0 + e] == k) { chain[j] = e; tab[h] = (int32_t)j; break; }
                h = (h + 1) & m;
            }
        }
        int64_t c = 0, o = (J->phase == 1) ? J->bcnt[b] : 0;
        for (int64_t i = a0; i < a1; i++) {
            uint32_t k = J->akey[i], h = hsh(k, J->pbits) & m;
            for (;;) {
                int32_t e = tab[h];
                if (e < 0) break;
                if (J->bkey[b0 + e] == k) {
                    for (int32_t x = e; x >= 0; x = chain[x]) {
                        if (J->phase == 1) { J->outA[o + c] = J->aidx[i]; J->outB[o + c] = J->bidx[b0 + x]; }
                        c++;
                    }
                    break;
                }
                h = (h + 1) & m;
            }
        }
        if (J->phase == 0) J->bcnt[b] = c;
    }
    free(tab); free(chain);
    return NULL;
}

int64_t kt_join(const uint32_t *ka, int64_t nA, const uint32_t *kb, int64_t nB,
                uint32_t *outA, uint32_t *outB, int64_t cap, int threads) {
    int nt = threads < 1 ? 1 : (threads > 256 ? 256 : threads);
    int pbits = 0; while (((int64_t)1 << pbits) * 2048 < nB && pbits < 24) pbits++;
    int64_t nb = (int64_t)1 << pbits;
    uint32_t *akey = malloc(4 * (nA + 1)), *aidx = malloc(4 * (nA + 1));
    uint32_t *bkey = malloc(4 * (nB + 1)), *bidx = malloc(4 * (nB + 1));
    int64_t *offA = malloc(8 * (nb + 1)), *offB = malloc(8 * (nb + 1)), *bcnt = malloc(8 * nb);
    partition(ka, nA, pbits, nt, akey, aidx, offA);
    partition(kb, nB, pbits, nt, bkey, bidx, offB);
    int64_t total = 0;
    pthread_t th[256]; join_t J[256];
    for (int ph = 0; ph < 2; ph++) {
        volatile int64_t next = 0;
        for (int t = 0; t < nt; t++) {
            J[t] = (join_t){akey, aidx, bkey, bidx, offA, offB, nb, pbits, ph, bcnt, outA, outB, &next};
            pthread_create(&th[t], NULL, join_worker, &J[t]);
        }
        for (int t = 0; t < nt; t++) pthread_join(th[t], NULL);
        if (ph == 0) {
            int64_t s = 0;
            for (int64_t b = 0; b < nb; b++) { int64_t c = bcnt[b]; bcnt[b] = s; s += c; }
            total = s;
            if (total > cap) break;
        }
    }
    free(akey); free(aidx); free(bkey); free(bidx); free(offA); free(offB); free(bcnt);
    return total;
}

typedef struct { const uint32_t *x, *y, *ia, *ib; int64_t lo, hi; uint32_t M; uint32_t *out; } mg_t;

static void *mg_worker(void *v) {
    mg_t *g = v;
    for (int64_t k = g->lo; k < g->hi; k++)
        g->out[k] = (uint32_t)((uint64_t)g->x[g->ia[k]] * g->y[g->ib[k]] % g->M);
    return NULL;
}

void kt_mulgather(const uint32_t *x, const uint32_t *y, const uint32_t *ia, const uint32_t *ib,
                  int64_t n, uint32_t M, uint32_t *out, int threads) {
    int nt = threads < 1 ? 1 : (threads > 256 ? 256 : threads);
    pthread_t th[256]; mg_t g[256];
    for (int t = 0; t < nt; t++) {
        g[t] = (mg_t){x, y, ia, ib, n * t / nt, n * (t + 1) / nt, M, out};
        pthread_create(&th[t], NULL, mg_worker, &g[t]);
    }
    for (int t = 0; t < nt; t++) pthread_join(th[t], NULL);
}
