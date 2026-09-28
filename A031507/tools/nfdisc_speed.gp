default(nbthreads,1);
t=getabstime(); L=nflist("S3",[10^7,10^7+10^6]); t1=getabstime()-t;
t=getabstime(); n=0; for(i=1,#L, my(d=nfdisc(L[i])); if(isfundamental(d), n++)); t2=getabstime()-t;
t=getabstime(); m=0; for(i=1,#L, my(p=poldisc(L[i])); if(isfundamental(p), m++)); t3=getabstime()-t;
printf("|disc| in [1e7,1.1e7]: %d fields; nflist %d ms; nfdisc loop %d ms (%d fundamental); poldisc-fundamental %d (%d ms)\n",#L,t1,t2,n,m,t3);
