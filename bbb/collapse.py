"""Exact necessary seventh-power phase systems for BBB triples at width 11.

This is a support/phase reduction, NOT a classification of original triples.
The seventh-root all-different lifting conditions are checked separately.
"""
from fractions import Fraction as Q
from itertools import combinations
from collections import Counter, deque
from pathlib import Path
import json, time
from torus_solver import solve_torus

BASE = (tuple(range(7)), (7,8), (9,10))

def partitions():
    out=[]
    for small in combinations(range(11),4):
        a,b,c,d=small
        big=tuple(j for j in range(11) if j not in small)
        for p,q in (((a,b),(c,d)),((a,c),(b,d)),((a,d),(b,c))):
            out.append((big,*sorted((p,q))))
    return sorted(out)

def act(P,g):
    return (tuple(sorted(g[x] for x in P[0])),
            *sorted(tuple(sorted(g[x] for x in b)) for b in P[1:]))

def orbits(parts):
    gens=[]
    for a,b in [(i,i+1) for i in range(6)]+[(7,8),(9,10)]:
        g=list(range(11)); g[a],g[b]=g[b],g[a]; gens.append(g)
    gens.append(list(range(7))+[9,10,7,8])
    left=set(parts); out=[]
    while left:
        p=min(left); seen={p}; todo=deque([p])
        while todo:
            r=todo.popleft()
            for g in gens:
                s=act(r,g)
                if s not in seen: seen.add(s); todo.append(s)
        assert seen <= left
        left-=seen; out.append((p,len(seen)))
    return out

def row(P,offset):
    c=[Q(0)]*11; B=[[0]*4 for _ in range(11)]
    for k,pair in enumerate(P[1:]):
        for j in pair: B[j][offset+k]=1
        c[pair[1]]=Q(1,2)
    return c,B

def equations(P,R):
    uc,U=row(BASE,0); vc,V=row(P,2)
    c=[a-b for a,b in zip(uc,vc)]
    B=[[a-b for a,b in zip(x,y)] for x,y in zip(U,V)]
    # All three 7-cycles have seventh power 1.  Their phases are compatible
    # because their three supports have a common column in every genuine BBB.
    A=[B[j] for j in R[0]]
    rhs=[-c[j] % 1 for j in R[0]]
    for a,b in R[1:]:
        A.append([y-x for x,y in zip(B[a],B[b])])
        rhs.append((Q(1,2)-c[b]+c[a]) % 1)
    return A,rhs

def main():
    t=time.perf_counter(); parts=partitions(); reps=orbits(parts)
    assert len(parts)==990 and sum(n for _,n in reps)==990
    report={'scope':'necessary collapsed BBB systems; lifting not yet imposed',
            'partition_count':len(parts),'support_orbits':[], 'records':[]}
    totals=Counter()
    for vi,(P,size) in enumerate(reps):
        count=Counter()
        for wi,R in enumerate(parts):
            count['systems']+=1
            if not set(BASE[0]) & set(P[0]) & set(R[0]):
                count['empty_triple_intersection']+=1; continue
            A,b=equations(P,R)
            sols=solve_torus(A,b,4)
            count['consistent_systems']+=bool(sols)
            for c,B in sols:
                dim=len(B[0]); count[f'dimension_{dim}']+=1
                for a,rhs in zip(A,b):
                    assert (sum(x*y for x,y in zip(a,c))-rhs)%1==0
                    assert all(sum(a[i]*B[i][k] for i in range(4))==0 for k in range(dim))
                report['records'].append({'vi':vi,'wi':wi,'c':list(map(str,c)),
                                          'B':B,'dimension':dim})
        report['support_orbits'].append({'representative':P,'size':size,'counts':dict(count)})
        totals.update(count)
        print(vi,P,size,dict(count),flush=True)
    report['totals']=dict(totals)
    report['constant_denominators']=sorted({Q(x).denominator for r in report['records'] for x in r['c']})
    report['elapsed_seconds']=time.perf_counter()-t
    Path(__file__).with_name('collapsed.json').write_text(json.dumps(report,indent=2)+'\n')
    print('TOTAL',dict(totals),'denominators',report['constant_denominators'],flush=True)

if __name__=='__main__': main()
