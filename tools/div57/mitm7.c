// Find subsets S of the pool with  f * prod(S) ≡ 1 (mod L), f a forced prime (7).
// Side A (na primes) in a hash table; side B (nb primes) streamed by Gray code over
// T threads. Usage: ./mitm7 forced na nb threads < pool file ("L k q1..qk").
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
typedef unsigned __int128 u128; typedef __int128 i128;
static uint64_t L; static int k; static uint64_t q[128], qi[128];
static int na, nb, T; static uint64_t forced;
static uint64_t *key; static uint64_t *amask; static uint64_t tmask;
static inline uint64_t mm(uint64_t a,uint64_t b){return (uint64_t)((u128)a*b%L);}
static uint64_t inv(uint64_t a){ i128 t=0,nt=1,r=L,nr=a%L; while(nr){ i128 qq=r/nr,tmp=t-qq*nt; t=nt; nt=tmp; tmp=r-qq*nr; r=nr; nr=tmp;} if(t<0)t+=L; return (uint64_t)t;}
static inline uint64_t h(uint64_t x){ x^=x>>31; x*=0x9E3779B97F4A7C15ULL; x^=x>>29; return x&tmask; }
static pthread_mutex_t mx=PTHREAD_MUTEX_INITIALIZER; static long nsol=0;
static void report(uint64_t ma, uint64_t mb){
  pthread_mutex_lock(&mx); nsol++;
  printf("SOL %llu",(unsigned long long)forced);
  for(int i=0;i<na;i++) if(ma>>i&1) printf(" %llu",(unsigned long long)q[i]);
  for(int i=0;i<nb;i++) if(mb>>i&1) printf(" %llu",(unsigned long long)q[na+i]);
  printf("\n"); fflush(stdout); pthread_mutex_unlock(&mx);
}
static void lookup(uint64_t y, uint64_t mb){
  uint64_t s=h(y); while(key[s]){ if(key[s]==y) report(amask[s],mb); s=(s+1)&tmask; }
}
typedef struct { int id; } Arg;
static void *work(void *vp){
  int id=((Arg*)vp)->id; int top=0; while((1<<top)<T) top++;
  int low=nb-top; // bits [0,low) gray-coded, bits [low,nb) fixed by id
  uint64_t y=inv(forced); uint64_t mb=0;
  for(int i=0;i<top;i++) if(id>>i&1){ y=mm(y,qi[na+low+i]); mb|=1ULL<<(low+i); }
  lookup(y,mb);
  uint64_t N=1ULL<<low;
  for(uint64_t g=1; g<N; g++){
    int j=__builtin_ctzll(g); uint64_t bit=1ULL<<j;
    if(mb&bit){ mb&=~bit; y=mm(y,q[na+j]); } else { mb|=bit; y=mm(y,qi[na+j]); }
    lookup(y,mb);
  }
  return 0;
}
int main(int argc,char**argv){
  forced=strtoull(argv[1],0,10); na=atoi(argv[2]); nb=atoi(argv[3]); T=atoi(argv[4]);
  unsigned long long LL; scanf("%llu %d",&LL,&k); L=LL; for(int i=0;i<k;i++){unsigned long long x; scanf("%llu",&x); q[i]=x;}
  for(int i=0;i<na+nb;i++) qi[i]=inv(q[i]);
  uint64_t NA=1ULL<<na; uint64_t tsz=1; while(tsz<2*NA) tsz<<=1; tmask=tsz-1;
  key=calloc(tsz,8); amask=malloc(8*tsz);
  // build side A by Gray code
  uint64_t x=1, m=0; 
  for(uint64_t g=0; g<NA; g++){
    if(g){ int j=__builtin_ctzll(g); uint64_t bit=1ULL<<j; if(m&bit){m&=~bit; x=mm(x,qi[j]);} else {m|=bit; x=mm(x,q[j]);} }
    uint64_t s=h(x); while(key[s]) s=(s+1)&tmask; key[s]=x; amask[s]=m;
  }
  fprintf(stderr,"side A built: 2^%d entries\n",na);
  pthread_t th[256]; Arg a[256];
  for(int i=0;i<T;i++){ a[i].id=i; pthread_create(&th[i],0,work,&a[i]); }
  for(int i=0;i<T;i++) pthread_join(th[i],0);
  fprintf(stderr,"done: solutions %ld\n",nsol);
}
