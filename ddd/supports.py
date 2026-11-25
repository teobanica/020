"""Chosen support partitions for D=5+2+2+2 at width 11."""
from itertools import combinations
from collections import deque
from pathlib import Path
import json
BASE=(tuple(range(5)),(5,6),(7,8),(9,10))

def matchings(rem):
    if not rem:
        yield ();return
    a=rem[0]
    for b in rem[1:]:
        for rest in matchings(tuple(j for j in rem[1:] if j!=b)):
            yield ((a,b),)+rest

def partitions():
    return sorted((big,)+small for big in combinations(range(11),5)
        for small in matchings(tuple(j for j in range(11) if j not in big)))

def act(P,g):
    return (tuple(sorted(g[j] for j in P[0])),
            *sorted(tuple(sorted(g[j] for j in block)) for block in P[1:]))

def orbits(parts):
    gens=[]
    for a,b in [(0,1),(1,2),(2,3),(3,4),(5,6),(7,8),(9,10)]:
        g=list(range(11));g[a],g[b]=g[b],g[a];gens.append(g)
    gens.append(list(range(5))+[7,8,5,6,9,10])
    gens.append(list(range(5))+[5,6,9,10,7,8])
    left=set(parts);out=[]
    while left:
        p=min(left);seen={p};todo=deque([p])
        while todo:
            r=todo.popleft()
            for g in gens:
                q=act(r,g)
                if q not in seen:seen.add(q);todo.append(q)
        assert seen<=left
        left-=seen;out.append((p,len(seen)))
    return out

def row(P,offset):
    B=[[0]*6 for _ in range(11)]
    for k,block in enumerate(P[1:]):
        for j in block:B[j][offset+k]=1
    return B

def offsets(P):
    out=[0]*11
    for block in P[1:]:out[block[1]]=1
    return out

def equations(Q,R):
    U,V=row(BASE,0),row(Q,3)
    W=[[a-b for a,b in zip(x,y)] for x,y in zip(U,V)]
    ep,eq,er=offsets(BASE),offsets(Q),offsets(R)
    h=[c-a+b for a,b,c in zip(ep,eq,er)]
    A=[];g=[]
    for block in R:
        ref=block[0]
        for j in block[1:]:
            A.append([a-b for a,b in zip(W[j],W[ref])]);g.append(h[j]-h[ref])
    return A,g

def main():
    root=Path(__file__).resolve().parent;parts=partitions();reps=orbits(parts)
    data={'partition_count':len(parts),'support_orbits':[{'representative':p,'size':n} for p,n in reps]}
    (root/'supports.json').write_text(json.dumps(data,indent=2)+'\n')
    with (root/'phase_input.txt').open('w') as f:
        f.write(str(len(parts)*len(reps))+'\n')
        for vi,(Q,_) in enumerate(reps):
            for wi,R in enumerate(parts):
                A,g=equations(Q,R);f.write(f'{vi} {wi}\n')
                for a,b in zip(A,g):f.write(' '.join(map(str,a+[b]))+'\n')
    print(len(parts),'partitions;',len(reps),'relative orbits;',len(parts)*len(reps),'systems',flush=True)
    print([n for _,n in reps],flush=True)

if __name__=='__main__':main()
