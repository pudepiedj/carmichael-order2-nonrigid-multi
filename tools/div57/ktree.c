// 4-list (Wagner-style) search for  f * prod(S) ≡ 1 (mod L1*L2),  gcd(L1,L2)=1, both < 2^64.
// Lists A,B,C,D of pool primes. Level 1: a*b ≡ t (mod L1) and c*d ≡ (f t)^-1 (mod L1).
// Level 2: match (ab)*(cd)*f ≡ 1 (mod L2).
// Input: "L1 L2 f t_seed nA nB nC nD q1 ... " ; argv[1] = threads.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
typedef unsigned __int128 u128; typedef __int128 i128;
static uint64_t L1,L2,f; static int nA,nB,nC,nD,T; static uint64_t q[256];
static inline uint64_t m1(uint64_t a,uint64_t b){return (uint64_t)((u128)a*b%L1);}
static inline uint64_t m2(uint64_t a,uint64_t b){return (uint64_t)((u128)a*b%L2);}
static uint64_t inv(uint64_t a,uint64_t m){ i128 t=0,nt=1,r=m,nr=a%m; while(nr){ i128 qq=r/nr,tmp=t-qq*nt; t=nt; nt=tmp; tmp=r-qq*nr; r=nr; nr=tmp;} if(t<0)t+=m; return (uint64_t)t;}
static inline uint64_t hh(uint64_t x,uint64_t mask){ x^=x>>31; x*=0x9E3779B97F4A7C15ULL; x^=x>>29; return x&mask; }
// hash for level 1: key (mod L1), payload (mod L2), mask
typedef struct { uint64_t key, pay; uint32_t mask; } H1;
static H1 *h1; static uint64_t h1m;
static void build_h1(int off,int n,int store_inverse){
  uint64_t N=1ULL<<n, sz=1; while(sz<2*N) sz<<=1; h1m=sz-1; h1=calloc(sz,sizeof(H1));
  uint64_t x1=1%L1,x2=1%L2,mask=0;
  for(uint64_t g=0; g<N; g++){
    if(g){int j=__builtin_ctzll(g); uint64_t b=1ULL<<j; uint64_t p=q[off+j];
      if(mask&b){mask&=~b; x1=m1(x1,inv(p,L1)); x2=m2(x2,inv(p,L2));} else {mask|=b; x1=m1(x1,p); x2=m2(x2,p);} }
    uint64_t s=hh(x1,h1m); while(h1[s].key) s=(s+1)&h1m;
    h1[s].key=x1; h1[s].pay= store_inverse? inv(x2,L2) : x2; h1[s].mask=(uint32_t)mask;
  }
}
// level-1 output
typedef struct { uint64_t v; uint32_t ma, mb; } P1;
typedef struct { P1 *buf; long n, cap; } Buf;
static Buf bufs[256];
static uint64_t tgt; static int phase; // phase 1: A x B -> list ; phase 2: C x D -> lookup
static P1 *list1; static long n1; static uint64_t *h2key; static uint32_t *h2idx; static uint64_t h2m;
static uint64_t qinv1[256], qinv2[256];
static pthread_mutex_t mx=PTHREAD_MUTEX_INITIALIZER; static long nsol=0;
typedef struct { int id; } Arg;
static void *stream(void *vp){
  int id=((Arg*)vp)->id; int off = phase==1? nA : nA+nB+nC; int n = phase==1? nB : nD;
  int top=0; while((1<<top)<T) top++; int low=n-top;
  // maintain prod of inverses mod L1 (for the key) and prod mod L2 (phase 1) / inverse mod L2 (phase 2)
  uint64_t y1=tgt, y2=1%L2, mask=0;
  for(int i=0;i<top;i++) if(id>>i&1){ int j=low+i; y1=m1(y1,qinv1[off+j]); y2= phase==1? m2(y2,q[off+j]) : m2(y2,qinv2[off+j]); mask|=1ULL<<j; }
  uint64_t N=1ULL<<low; Buf *bf=&bufs[id];
  for(uint64_t g=0; g<N; g++){
    if(g){int j=__builtin_ctzll(g); uint64_t b=1ULL<<j; int k=off+j;
      if(mask&b){mask&=~b; y1=m1(y1,q[k]); y2= phase==1? m2(y2,qinv2[k]) : m2(y2,q[k]);}
      else {mask|=b; y1=m1(y1,qinv1[k]); y2= phase==1? m2(y2,q[k]) : m2(y2,qinv2[k]);} }
    uint64_t s=hh(y1,h1m);
    while(h1[s].key){ if(h1[s].key==y1){
        if(phase==1){ if(bf->n==bf->cap){bf->cap=bf->cap?2*bf->cap:1<<20; bf->buf=realloc(bf->buf,bf->cap*sizeof(P1));}
          bf->buf[bf->n].v=m2(h1[s].pay,y2); bf->buf[bf->n].ma=h1[s].mask; bf->buf[bf->n].mb=(uint32_t)mask; bf->n++; }
        else { // need ab ≡ f^-1 * c^-1 * d^-1 (mod L2)
          uint64_t key=m2(m2(inv(f,L2),h1[s].pay),y2);
          uint64_t s2=hh(key,h2m);
          while(h2key[s2]){ if(h2key[s2]==key){ P1 *p=&list1[h2idx[s2]];
              pthread_mutex_lock(&mx); nsol++; printf("SOL %llu",(unsigned long long)f);
              for(int i=0;i<nA;i++) if(p->ma>>i&1) printf(" %llu",(unsigned long long)q[i]);
              for(int i=0;i<nB;i++) if(p->mb>>i&1) printf(" %llu",(unsigned long long)q[nA+i]);
              for(int i=0;i<nC;i++) if(h1[s].mask>>i&1) printf(" %llu",(unsigned long long)q[nA+nB+i]);
              for(int i=0;i<nD;i++) if(mask>>i&1) printf(" %llu",(unsigned long long)q[nA+nB+nC+i]);
              printf("\n"); fflush(stdout); pthread_mutex_unlock(&mx); }
            s2=(s2+1)&h2m; }
        } }
      s=(s+1)&h1m; }
  }
  return 0;
}
static void run_stream(){ pthread_t th[256]; Arg a[256]; for(int i=0;i<T;i++){a[i].id=i; pthread_create(&th[i],0,stream,&a[i]);} for(int i=0;i<T;i++) pthread_join(th[i],0); }
int main(int argc,char**argv){
  T=atoi(argv[1]);
  unsigned long long a1,a2,ff,seed; scanf("%llu %llu %llu %llu %d %d %d %d",&a1,&a2,&ff,&seed,&nA,&nB,&nC,&nD); L1=a1; L2=a2; f=ff;
  int k=nA+nB+nC+nD; for(int i=0;i<k;i++){unsigned long long x; scanf("%llu",&x); q[i]=x; qinv1[i]=inv(x,L1); qinv2[i]=inv(x,L2);}
  // target t: product of the pool primes selected by seed bits, mod L1 (a unit)
  uint64_t t=1%L1; for(int i=0;i<64 && i<k;i++) if(seed>>i&1) t=m1(t,q[i]);
  // phase 1
  build_h1(0,nA,0); fprintf(stderr,"hash A built\n");
  phase=1; tgt=t; run_stream();
  n1=0; for(int i=0;i<T;i++) n1+=bufs[i].n; list1=malloc(n1*sizeof(P1)); long c=0;
  for(int i=0;i<T;i++){ memcpy(list1+c,bufs[i].buf,bufs[i].n*sizeof(P1)); c+=bufs[i].n; free(bufs[i].buf); bufs[i].buf=0; bufs[i].n=bufs[i].cap=0; }
  free(h1); fprintf(stderr,"level-1 list AB: %ld pairs\n",n1);
  uint64_t sz=1; while(sz<2*(uint64_t)n1) sz<<=1; h2m=sz-1; h2key=calloc(sz,8); h2idx=malloc(4*sz);
  for(long i=0;i<n1;i++){ uint64_t s=hh(list1[i].v,h2m); while(h2key[s]) s=(s+1)&h2m; h2key[s]=list1[i].v; h2idx[s]=(uint32_t)i; }
  fprintf(stderr,"level-2 hash built\n");
  // phase 2: c*d ≡ (f t)^-1 (mod L1)  =>  key c ≡ (f t)^-1 * d^-1
  build_h1(nA+nB,nC,1); fprintf(stderr,"hash C built\n");
  phase=2; tgt=inv(m1(f%L1,t),L1); run_stream();
  fprintf(stderr,"done: solutions %ld\n",nsol);
}
