"""Nonzero polynomial witnesses for the fourteen seed-degenerate presentations."""
from itertools import permutations,product
from pathlib import Path
import json
P=421;ROOT=2
PW=[pow(ROOT,i,P) for i in range(420)]
def det(A):
    n=len(A);ans=0
    for s in permutations(range(n)):
        inv=sum(s[i]>s[j] for i in range(n) for j in range(i+1,n));v=(-1)**inv
        for i,j in enumerate(s):v=v*A[i][j]%P
        ans+=v
    return ans%P
H=json.loads(Path('seed.json').read_text())['exponents'];out=[]
for idx,line in enumerate(Path('unresolved.txt').read_text().splitlines()):
    t=line.split();b=list(map(int,t[1:12]));e=list(map(int,t[12:]))
    for shifts in product((0,5,10,15),repeat=4):
        A=[];B=[]
        for r in range(4):
            ex=[(e[j]-21*(H[r+1][j]+(shifts[r] if j>=5 else 0)))%420 for j in range(11)]
            A.append([sum(PW[ex[j]] for j in range(11) if b[j]==k)%P for k in range(5)])
            B.append([sum(PW[-ex[j]%420] for j in range(11) if b[j]==k)%P for k in range(5)])
        norms=[det([[row[j] for j in range(5) if j!=omit] for row in A])*det([[row[j] for j in range(5) if j!=omit] for row in B])%P for omit in range(5)]
        if len(set(norms))>1:
            out.append({'index':idx,'outer_phase_exponents_mod20':shifts,'cofactor_norms_mod421':norms});break
    else:raise AssertionError(('no generic witness',idx))
Path('generic_witnesses.json').write_text(json.dumps(out,indent=2)+'\n')
print('All 14 have nonzero Laurent-polynomial obstructions. Witness phase tuples:',sorted(set(tuple(x['outer_phase_exponents_mod20']) for x in out)))
