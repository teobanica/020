"""Support partitions and exact phase equations for CCC triples, width 11."""
from itertools import combinations
from collections import deque
from pathlib import Path
import json

BASE=(tuple(range(5)),(5,6,7),(8,9,10))

def partitions():
    out=[]
    for big in combinations(range(11),5):
        rem=tuple(x for x in range(11) if x not in big)
        for extra in combinations(rem[1:],2):
            a=(rem[0],)+extra
            b=tuple(x for x in rem if x not in a)
            out.append((big,a,b))
    return sorted(out)

def act(P,g):
    return (tuple(sorted(g[x] for x in P[0])),
            *sorted(tuple(sorted(g[x] for x in block)) for block in P[1:]))

def orbits(parts):
    gens=[]
    for a,b in [(0,1),(1,2),(2,3),(3,4),(5,6),(6,7),(8,9),(9,10)]:
        g=list(range(11));g[a],g[b]=g[b],g[a];gens.append(g)
    gens.append(list(range(5))+[8,9,10,5,6,7])
    left=set(parts);out=[]
    while left:
        p=min(left);seen={p};todo=deque([p])
        while todo:
            for g in gens:
                q=act(todo[0],g)
                if q not in seen:seen.add(q);todo.append(q)
            todo.popleft()
        assert seen<=left
        left-=seen;out.append((p,len(seen)))
    return out

def row(P,offset):
    B=[[0]*4 for _ in range(11)]
    for k,block in enumerate(P[1:]):
        for j in block:B[j][offset+k]=1
    return B

def equations(Q,R):
    U,V=row(BASE,0),row(Q,2)
    W=[[a-b for a,b in zip(x,y)] for x,y in zip(U,V)]
    return [[a-b for a,b in zip(W[j],W[block[0]])]
            for block in R for j in block[1:]]

def main():
    root=Path(__file__).resolve().parent
    parts=partitions();reps=orbits(parts)
    data={'partition_count':len(parts),'support_orbits':[{'representative':p,'size':n} for p,n in reps]}
    (root/'supports.json').write_text(json.dumps(data,indent=2)+'\n')
    with (root/'phase_input.txt').open('w') as f:
        f.write(str(len(reps)*len(parts))+'\n')
        for vi,(Q,_) in enumerate(reps):
            for wi,R in enumerate(parts):
                A=equations(Q,R)
                f.write(f'{vi} {wi} {len(A)}\n')
                for a in A:f.write(' '.join(map(str,a))+'\n')
    print(len(parts),'partitions;',len(reps),'relative orbits;',len(parts)*len(reps),'systems',flush=True)
    print([n for _,n in reps],flush=True)

if __name__=='__main__':main()
