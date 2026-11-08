a(n)={my(h=(n+1)\2, id=matrix(h,h), k=0, c=0, o=i->min(i,n+1-i), M); for(b=1,h, for(i=1,b, id[i,b]=k++)); M=matrix(n,n,i,j,my(x=o(i),y=o(j)); id[min(x,y),max(x,y)]); for(m=0,2^k-1, my(A=matrix(n,n,i,j,bittest(m,M[i,j]-1))); if(prod(i=1,n,prod(j=1,n, !A[i,j] || sum(u=max(i-1,1),min(i+1,n), sum(v=max(j-1,1),min(j+1,n), A[u,v]))>1)), c++)); c};
for(n=1,10, print(n," ",a(n)));
