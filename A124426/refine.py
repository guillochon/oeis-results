from collections import defaultdict
from math import comb
from sympy import factorint
def refined(N):
    # state (s,f,e,t,k) -> weight ; k = #open dolls
    st={(0,0,0,0,0):1}; res=[]
    def place(states, kind):
        ns=defaultdict(int)
        for (s,f,e,t,k),w in states.items():
            opts=[(1,0,0,0,0,False),(s,-1,0,0,0,False)]
            opts+= [(f,0,-1,0,0,False),(e,0,0,-1,0,True)] if kind in 'DB' else [(t,0,0,0,-1,True)]
            for m,ds,df,de,dt,enc in opts:
                if not m: continue
                s2,f2,e2,t2=s+ds,f+df,e+de,t+dt
                if kind=='D': e2+=1; s2+=(not enc)
                if kind=='B':
                    if enc: e2+=1
                    else: f2+=1
                if kind=='T': t2+=1; s2+=(not enc)
                ns[(s2,f2,e2,t2,k+(kind=='T'))]+=w*m
        return ns
    for n in range(1,N+1):
        a=place(st,'D'); b=place(place(st,'B'),'T')
        st=defaultdict(int)
        for d in (a,b):
            for x,v in d.items(): st[x]+=v
        row=defaultdict(int)
        for x,v in st.items(): row[x[4]]+=v
        res.append([row[k] for k in range(n+1)])
    return res
for n,row in enumerate(refined(7),1):
    print(n,row,[r//comb(n,k) if r%comb(n,k)==0 else None for k,r in enumerate(row)])
