/* Random D4-symmetric n x n designs with no pinhole (open cell, all 8 neighbors filled; outside = filled),
   density-matched: Metropolis over symmetry orbits with a field h tuned during burn-in so the mean fill
   approaches the target.  Usage: null_sampler n target seed nsamples  -> one line of n*n 0/1 per sample. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static int n, A[80][80], orb[80][80], K, ocell[2000][8][2], ocnt[2000];
static unsigned long long s;
static double rnd(void){ s^=s<<13; s^=s>>7; s^=s<<17; return (s>>11)*(1.0/9007199254740992.0); }
static int at(int i,int j){ return (i<0||j<0||i>=n||j>=n)?1:A[i][j]; }
static int pin(int i,int j){ if(i<0||j<0||i>=n||j>=n||A[i][j]) return 0;
  for(int a=-1;a<=1;a++)for(int b=-1;b<=1;b++) if((a||b)&&!at(i+a,j+b)) return 0; return 1; }
static int okaround(int k){ for(int c=0;c<ocnt[k];c++){int i=ocell[k][c][0],j=ocell[k][c][1];
  for(int a=-1;a<=1;a++)for(int b=-1;b<=1;b++) if(pin(i+a,j+b)) return 0;} return 1; }
int main(int argc,char**argv){
  n=atoi(argv[1]); double target=atof(argv[2]); s=0x9E3779B97F4A7C15ULL*(atoi(argv[3])+1); int ns=atoi(argv[4]);
  int id[80][80]; for(int i=0;i<80;i++)for(int j=0;j<80;j++) id[i][j]=-1; K=0;
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){int x=i<n-1-i?i:n-1-i,y=j<n-1-j?j:n-1-j,a=x<y?x:y,b=x<y?y:x;
    if(id[a][b]<0) id[a][b]=K++; int k=id[a][b]; orb[i][j]=k; ocell[k][ocnt[k]][0]=i; ocell[k][ocnt[k]][1]=j; ocnt[k]++; A[i][j]=0;}
  /* start all filled (valid) */
  for(int i=0;i<n;i++)for(int j=0;j<n;j++) A[i][j]=1;
  double h=0; long burn=4000L*K, gap=200L*K;
  for(int smp=0; smp<ns; smp++){
    long steps = smp==0?burn:gap;
    for(long t=0;t<steps;t++){
      int k=(int)(rnd()*K); int c0=ocell[k][0][0], c1=ocell[k][0][1]; int v=A[c0][c1], nv=1-v;
      double dE = (nv-v)*ocnt[k]*h; if(dE<0 && rnd()>exp(dE)) continue;
      for(int c=0;c<ocnt[k];c++) A[ocell[k][c][0]][ocell[k][c][1]]=nv;
      if(!okaround(k)) for(int c=0;c<ocnt[k];c++) A[ocell[k][c][0]][ocell[k][c][1]]=v;
      if(smp==0 && t%(20*K)==0){ int f=0; for(int i=0;i<n;i++)for(int j=0;j<n;j++) f+=A[i][j];
        h += 0.5*(target - (double)f/(n*n)) * (t<burn/2?1.0:0.3); }
    }
    for(int i=0;i<n;i++)for(int j=0;j<n;j++) putchar('0'+A[i][j]); putchar('\n');
  }
  fprintf(stderr,"h=%g\n",h); return 0; }
