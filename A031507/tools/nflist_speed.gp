default(nbthreads,1);
default(parisizemax, 4*10^9);
for(e=4,7, my(t=getabstime(), L=nflist("S3",[-10^e,-1])); printf("neg to 1e%d: %d fields  %d ms\n", e, #L, getabstime()-t));
for(e=4,7, my(t=getabstime(), L=nflist("S3",[1,10^e])); printf("pos to 1e%d: %d fields  %d ms\n", e, #L, getabstime()-t));
