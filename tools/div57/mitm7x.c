// mitm7x: mitm7.c for moduli D beyond 2^64.
//
// Finds every subset S of the pool with  f * prod(S) ≡ 1 (mod D),  f a forced prime.
// D is given as two coprime factors D1, D2 < 2^63, and every residue is kept as the pair
// (x mod D1, x mod D2). By the CRT the pair determines x mod D, so a match of pairs is a
// match mod D. Products stay in 64 bits (via 128-bit multiplication) for any D < 2^126.
//
// Side A (the first na primes) goes into a hash table: 2^(na+1) slots of 24 bytes.
// Side B (the next nb primes) is streamed by Gray code over T threads.
//
// Usage: ./mitm7x forced na nb threads < pool file ("D1 D2 k q1 .. qk")
// Checks: na + nb <= k (warns if < k: primes beyond na + nb are NOT searched);
// threads a power of two, at most 2^nb. Prints timings to stderr.
// Output: "SOL f q ..." lines on stdout, one per solution.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <time.h>
typedef unsigned __int128 u128; typedef __int128 i128;
static uint64_t D1, D2; static int k; static uint64_t q[128], qi1[128], qi2[128], q1[128], q2[128];
static int na, nb, T; static uint64_t forced;
static uint64_t *key1, *key2, *amask; static uint64_t tmask;
static inline uint64_t mm(uint64_t a, uint64_t b, uint64_t m){ return (uint64_t)((u128)a*b % m); }
static uint64_t inv(uint64_t a, uint64_t m){
  i128 t=0, nt=1, r=m, nr=a%m;
  while(nr){ i128 qq=r/nr, tmp=t-qq*nt; t=nt; nt=tmp; tmp=r-qq*nr; r=nr; nr=tmp; }
  if(r!=1){ fprintf(stderr,"%llu is not invertible mod %llu\n",(unsigned long long)a,(unsigned long long)m); exit(1); }
  if(t<0) t+=m; return (uint64_t)t;
}
static uint64_t gcd(uint64_t a, uint64_t b){ while(b){ uint64_t t=a%b; a=b; b=t; } return a; }
static inline uint64_t h(uint64_t x, uint64_t y){
  x ^= y*0xD6E8FEB86659FD93ULL; x ^= x>>31; x *= 0x9E3779B97F4A7C15ULL; x ^= x>>29; return x & tmask;
}
static double now(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); return ts.tv_sec+1e-9*ts.tv_nsec; }
static pthread_mutex_t mx=PTHREAD_MUTEX_INITIALIZER; static long nsol=0;
static void report(uint64_t ma, uint64_t mb){
  pthread_mutex_lock(&mx); nsol++;
  printf("SOL %llu",(unsigned long long)forced);
  for(int i=0;i<na;i++) if(ma>>i&1) printf(" %llu",(unsigned long long)q[i]);
  for(int i=0;i<nb;i++) if(mb>>i&1) printf(" %llu",(unsigned long long)q[na+i]);
  printf("\n"); fflush(stdout); pthread_mutex_unlock(&mx);
}
static inline void lookup(uint64_t y1, uint64_t y2, uint64_t mb){
  uint64_t s=h(y1,y2);
  while(key1[s]){ if(key1[s]==y1 && key2[s]==y2) report(amask[s],mb); s=(s+1)&tmask; }
}
typedef struct { int id; int top; } Arg;
static void *work(void *vp){
  int id=((Arg*)vp)->id, top=((Arg*)vp)->top;
  int low=nb-top;                         // bits [0,low) Gray-coded, bits [low,nb) fixed by id
  uint64_t y1=inv(forced%D1,D1), y2=inv(forced%D2,D2), mb=0;
  for(int i=0;i<top;i++) if(id>>i&1){ y1=mm(y1,qi1[na+low+i],D1); y2=mm(y2,qi2[na+low+i],D2); mb|=1ULL<<(low+i); }
  lookup(y1,y2,mb);
  uint64_t N=1ULL<<low;
  for(uint64_t g=1; g<N; g++){
    int j=__builtin_ctzll(g); uint64_t bit=1ULL<<j;
    if(mb&bit){ mb&=~bit; y1=mm(y1,q1[na+j],D1); y2=mm(y2,q2[na+j],D2); }
    else      { mb|=bit;  y1=mm(y1,qi1[na+j],D1); y2=mm(y2,qi2[na+j],D2); }
    lookup(y1,y2,mb);
  }
  return 0;
}
int main(int argc, char **argv){
  if(argc!=5){ fprintf(stderr,"usage: %s forced na nb threads < pool (\"D1 D2 k q1..qk\")\n",argv[0]); return 1; }
  forced=strtoull(argv[1],0,10); na=atoi(argv[2]); nb=atoi(argv[3]); T=atoi(argv[4]);
  unsigned long long a1,a2;
  if(scanf("%llu %llu %d",&a1,&a2,&k)!=3){ fprintf(stderr,"bad pool header\n"); return 1; }
  D1=a1; D2=a2;
  if(D1<2 || D2<2 || D1>>63 || D2>>63 || gcd(D1,D2)!=1){ fprintf(stderr,"need coprime D1, D2 in [2, 2^63)\n"); return 1; }
  if(k>128){ fprintf(stderr,"k = %d > 128\n",k); return 1; }
  for(int i=0;i<k;i++){ unsigned long long x; if(scanf("%llu",&x)!=1){ fprintf(stderr,"pool has fewer than k primes\n"); return 1; } q[i]=x; }
  if(na<0 || nb<1 || na+nb>k || na>40 || nb>62){ fprintf(stderr,"need 0 <= na <= 40, 1 <= nb, na + nb <= k = %d\n",k); return 1; }
  if(na+nb<k) fprintf(stderr,"WARNING: na + nb = %d < k = %d: the last %d primes are not searched\n",na+nb,k,k-na-nb);
  int top=0; while((1<<top)<T) top++;
  if((1<<top)!=T || top>nb){ fprintf(stderr,"threads must be a power of two, at most 2^nb\n"); return 1; }
  for(int i=0;i<na+nb;i++){
    q1[i]=q[i]%D1; q2[i]=q[i]%D2; qi1[i]=inv(q1[i],D1); qi2[i]=inv(q2[i],D2);
  }
  double t0=now();
  uint64_t NA=1ULL<<na, tsz=1; while(tsz<2*NA) tsz<<=1; tmask=tsz-1;
  key1=calloc(tsz,8); key2=malloc(8*tsz); amask=malloc(8*tsz);
  if(!key1||!key2||!amask){ fprintf(stderr,"out of memory for 2^%d slots\n",__builtin_ctzll(tsz)); return 1; }
  uint64_t x1=1%D1, x2=1%D2, m=0;          // side A by Gray code
  for(uint64_t g=0; g<NA; g++){
    if(g){ int j=__builtin_ctzll(g); uint64_t bit=1ULL<<j;
      if(m&bit){ m&=~bit; x1=mm(x1,qi1[j],D1); x2=mm(x2,qi2[j],D2); }
      else     { m|=bit;  x1=mm(x1,q1[j],D1);  x2=mm(x2,q2[j],D2); } }
    uint64_t s=h(x1,x2); while(key1[s]) s=(s+1)&tmask; key1[s]=x1; key2[s]=x2; amask[s]=m;
  }
  double t1=now();
  fprintf(stderr,"side A built: 2^%d entries, %.1f s\n",na,t1-t0);
  pthread_t th[256]; Arg ar[256];
  for(int i=0;i<T;i++){ ar[i].id=i; ar[i].top=top; pthread_create(&th[i],0,work,&ar[i]); }
  for(int i=0;i<T;i++) pthread_join(th[i],0);
  double t2=now();
  fprintf(stderr,"side B streamed: 2^%d subsets, %.1f s\ndone: solutions %ld, total %.1f s\n",nb,t2-t1,nsol,t2-t0);
  return 0;
}
