"""Nesting only (no stacking): brute-force enumeration of every arrangement. Counts A124426(n) = Bell(n)*Bell(n+1).
unit = (kind, i, contents): kind D (closed doll), B (bottom half), T (top half); contents holds at most one unit."""
import sys
from functools import lru_cache
from itertools import combinations
# pieces: ('B',i) / ('T',i). unit = (kind, i, contents) kind in D,B,T ; contents = tuple(sorted units)
CAP = {'table':None,'D':1,'B':1,'T':1}  # None = unlimited
def subsets(s):
    s=sorted(s)
    for r in range(len(s)+1):
        for c in combinations(s,r): yield frozenset(c)
@lru_cache(None)
def level(S, mode):
    """all tuples of units using exactly pieces S, placed in container 'mode'"""
    if not S: return [()]
    cap = CAP[mode]
    i = max(p[1] for p in S)
    out=set()
    mine=[p for p in S if p[1]==i]
    smaller=frozenset(p for p in S if p[1]<i)
    def opts_outer():
        # yields (unit, used pieces)
        if len(mine)==2:
            for X in subsets(smaller):
                for c in level(X,'D'): yield (('D',i,c), X|frozenset(mine))
        for p in mine:
            for X in subsets(smaller):
                for c in level(X,p[0]): yield ((p[0],i,c), X|{p})
    def allowed(u):
        if mode=='T': return u[0]=='T'
        if mode in ('B','D'): return u[0]!='T'
        return True
    # choose outer unit(s) for pieces of index i, possibly both halves as separate units
    results=[]
    def rec(units, rem, todo):
        if not todo:
            for r in level(rem, mode): results.append(tuple(sorted(units+list(r))))
            return
        p=todo[0]
        if any(p in u_used for _,u_used in []): pass
        for X in subsets(rem & smaller):
            for c in level(X,p[0]): rec(units+[(p[0],i,c)], rem-X-{p}, todo[1:])
    if len(mine)==2:
        for X in subsets(smaller):
            for c in level(X,'D'): 
                for r in level(smaller-X, mode): results.append(tuple(sorted([('D',i,c)]+list(r))))
    rec([], S, mine)
    for t in results:
        if not all(allowed(u) for u in t): continue
        if mode=='B':
            tops=[u for u in t if u[0]=='T']
            if any(('B',u[1]) not in [(v[0],v[1]) for v in t] for u in tops): continue
        if cap is not None and len(t)>cap: continue
        out.add(t)
    return sorted(out)
def a(n):
    S=frozenset([('B',k) for k in range(1,n+1)]+[('T',k) for k in range(1,n+1)])
    return level(S,'table')
if __name__=='__main__':
    for n in range(1,int(sys.argv[1])+1): print(n,len(a(n)),flush=True)
