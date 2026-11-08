// Count n x n 0/1 designs with no pinhole: no 0 cell whose 8 neighbors are all 1 (outside the grid = 1, the frame).
// mode 0: designs symmetric under all 8 symmetries of the square; mode 1: all designs. Brute force.
#include <stdio.h>
#include <stdlib.h>
static int n, N, cellOf[4096];
static int ok(const unsigned char *g){
    for(int r=0;r<n;r++)for(int q=0;q<n;q++){ if(g[r*n+q])continue; int all=1;
        for(int dr=-1;dr<=1&&all;dr++)for(int dq=-1;dq<=1;dq++){ if(!dr&&!dq)continue; int rr=r+dr,qq=q+dq;
            if(rr>=0&&qq>=0&&rr<n&&qq<n&&!g[rr*n+qq]){all=0;break;}}
        if(all)return 0;}
    return 1;
}
int main(int argc,char**argv){
    int mode=atoi(argv[1]); n=atoi(argv[2]); N=n*n; unsigned char g[4096]; long long c=0;
    if(mode==0){ int h=(n+1)/2,k=0,id[64][64];
        for(int a=0;a<h;a++)for(int b=a;b<h;b++)id[a][b]=k++;
        for(int i=0;i<n;i++)for(int j=0;j<n;j++){int a=i<n-1-i?i:n-1-i,b=j<n-1-j?j:n-1-j; if(a>b){int t=a;a=b;b=t;} cellOf[i*n+j]=id[a][b];}
        for(long long m=0;m<(1LL<<k);m++){for(int i=0;i<N;i++)g[i]=(m>>cellOf[i])&1; c+=ok(g);}
    } else for(long long m=0;m<(1LL<<N);m++){for(int i=0;i<N;i++)g[i]=(m>>i)&1; c+=ok(g);}
    printf("%d %d %lld\n",mode,n,c); return 0;
}
