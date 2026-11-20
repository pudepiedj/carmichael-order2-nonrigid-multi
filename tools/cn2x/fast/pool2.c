/* pool2.c -- Howe's pool P(2, L) = { primes q, q not dividing L, q^2 - 1 | L }, without storing divisors.
 *
 *   pool2 THREADS p1 e1 p2 e2 ...        (the factorisation of L; prints the pool primes, one per line, unsorted)
 *
 * Walks every even divisor d = q - 1 of L below 2^62 depth-first (no list), and keeps q = d + 1 when
 * d + 2 divides L with v_p(d) + v_p(d + 2) <= e_p for every p, q is prime (deterministic Miller-Rabin for
 * 64 bits) and q does not divide L.  Threads take the work in slices of the first two prime powers.
 * Matches cn2_nonrigid_scout.pool() (checked on the pool-689 generator); used by cn2_rich_multi_scout.py.
 *
 * build: clang -O3 -mcpu=native -o pool2 pool2.c -lpthread
 */
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static int np; static uint64_t P[64]; static int E[64];
static const uint64_t LIM = (uint64_t)1 << 62;

static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) { return (unsigned __int128)a * b % m; }
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m) {
    uint64_t r = 1; a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}
static int isprime64(uint64_t n) {
    if (n < 2) return 0;
    static const uint64_t sp[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int i = 0; i < 12; i++) { if (n % sp[i] == 0) return n == sp[i]; }
    uint64_t d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        uint64_t x = powmod(sp[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int ok = 0;
        for (int r = 1; r < s; r++) { x = mulmod(x, x, n); if (x == n - 1) { ok = 1; break; } }
        if (!ok) return 0;
    }
    return 1;
}

typedef struct { int tid, nt; uint64_t *out; int64_t n, cap; } job_t;

static void test(uint64_t d, const int *v, job_t *J) {
    uint64_t m = d + 2;
    for (int i = 0; i < np; i++) {
        int w = 0;
        while (m % P[i] == 0) { m /= P[i]; w++; }
        if (v[i] + w > E[i]) return;
    }
    if (m != 1) return;
    uint64_t q = d + 1;
    for (int i = 0; i < np; i++) if (q == P[i]) return;
    if (!isprime64(q)) return;
    if (J->n == J->cap) { J->cap = J->cap ? 2 * J->cap : 1024; J->out = realloc(J->out, 8 * J->cap); }
    J->out[J->n++] = q;
}

static void dfs(int i, uint64_t d, int *v, job_t *J) {
    if (i == np) { if (d % 2 == 0) test(d, v, J); return; }
    uint64_t x = d;
    for (int k = 0; k <= E[i]; k++) {
        v[i] = k;
        dfs(i + 1, x, v, J);
        if (k < E[i]) { if (x > LIM / P[i]) break; x *= P[i]; }
    }
    v[i] = 0;
}

static void *worker(void *a) {
    job_t *J = a;
    int v[64] = {0};
    /* slice: the exponents of the first two primes, (k0, k1), task id k0*(E1+1)+k1, plus the third prime if present */
    int n0 = E[0] + 1, n1 = np > 1 ? E[1] + 1 : 1, n2 = np > 2 ? E[2] + 1 : 1;
    for (int t = J->tid; t < n0 * n1 * n2; t += J->nt) {
        int k0 = t / (n1 * n2), k1 = (t / n2) % n1, k2 = t % n2;
        uint64_t d = 1; int ok = 1;
        int ks[3] = {k0, k1, k2};
        for (int j = 0; j < 3 && j < np; j++)
            for (int r = 0; r < ks[j]; r++) { if (d > LIM / P[j]) { ok = 0; break; } d *= P[j]; }
        if (!ok) continue;
        for (int j = 0; j < 3 && j < np; j++) v[j] = ks[j];
        dfs(np < 3 ? np : 3, d, v, J);
    }
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 4 || (argc - 2) % 2) { fprintf(stderr, "usage: pool2 THREADS p1 e1 p2 e2 ...\n"); return 1; }
    int nt = atoi(argv[1]);
    np = (argc - 2) / 2;
    for (int i = 0; i < np; i++) { P[i] = strtoull(argv[2 + 2 * i], 0, 10); E[i] = atoi(argv[3 + 2 * i]); }
    pthread_t th[256]; job_t J[256];
    for (int t = 0; t < nt; t++) { J[t] = (job_t){t, nt, NULL, 0, 0}; pthread_create(&th[t], NULL, worker, &J[t]); }
    for (int t = 0; t < nt; t++) {
        pthread_join(th[t], NULL);
        for (int64_t k = 0; k < J[t].n; k++) printf("%llu\n", (unsigned long long)J[t].out[k]);
    }
    return 0;
}
