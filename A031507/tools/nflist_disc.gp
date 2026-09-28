default(nbthreads,1);
L=nflist("S3",[1,10^5]);
bad=0; for(i=1,#L, if(poldisc(L[i])!=nfdisc(L[i]), bad++));
print("fields ",#L," with poldisc != nfdisc: ",bad);
print(L[1..5]);
\ h3 check against quadclassunit for fundamental D with |D| <= 10^5
cnt=Map(); for(i=1,#L, my(d=nfdisc(L[i])); if(isfundamental(d), mapput(cnt,d,if(mapisdefined(cnt,d,&c),c+1,1))));
h3(D) = my(c=quadclassunit(D).cyc); sum(i=1,#c,c[i]%3==0);
mis=0; tot=0; forstep(D=-10^5,10^5,1, if(D!=0 && D!=1 && isfundamental(D), tot++; my(n=0); mapisdefined(cnt,D,&n); if(2*n+1 != 3^h3(D), mis++; if(mis<5,print("mismatch D=",D," n=",n," h3=",h3(D))))));
print("fundamental D checked: ",tot,", mismatches: ",mis);
