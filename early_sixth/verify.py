"""Check saved exact certificates without repeating their RREF solver."""
import json
from pathlib import Path
from collections import Counter
from exact_remaining import K,POW,PHI,divmodp,Q

def cyclotomic(n):
    f=[Q(-1)]+[Q(0)]*(n-1)+[Q(1)]
    for d in range(1,n):
        if n%d==0:
            f,r=divmodp(f,cyclotomic(d));assert not r
    return f
assert cyclotomic(60)==PHI
H=json.loads(Path('seed.json').read_text())['exponents']
assert H==[
 [0,0,0,0,0,0,0,0,0,0,0],
 [0,4,8,12,16,0,10,15,5,5,15],
 [0,8,16,4,12,0,15,10,15,5,5],
 [0,12,4,16,8,0,5,15,10,15,5],
 [0,16,12,8,4,0,5,5,15,10,15]]
assert all(sum((POW[3*(H[r][j]-H[s][j])%60] for j in range(11)),K())==0 for r in range(5) for s in range(r))
for cert in json.loads(Path('seed.json').read_text())['quotient_certificates']:
    r,s=cert['rows'];cyc=[cert['pentagon']]+cert['pairs']
    assert sorted(j for block in cyc for j in block)==list(range(11))
    for block in cyc:
        ex=[(H[r][j]-H[s][j])%20 for j in block];origin=ex[0]
        assert sorted((e-origin)%20 for e in ex)==list(range(0,20,20//len(block)))

def modular_rank(a,p=421):
    a=[row[:] for row in a];r=0
    for c in range(len(a[0])):
        i=next((i for i in range(r,len(a)) if a[i][c]%p),None)
        if i is None:continue
        a[r],a[i]=a[i],a[r];v=pow(a[r][c],-1,p);a[r]=[x*v%p for x in a[r]]
        for i in range(r+1,len(a)):
            v=a[i][c];a[i]=[(x-v*y)%p for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r

lines=Path('unresolved.txt').read_text().splitlines();records=json.loads(Path('exact_remaining.json').read_text())
assert len(lines)==len(records)==14
for idx,(line,rec) in enumerate(zip(lines,records)):
    assert rec['index']==idx
    w=line.split();b=list(map(int,w[1:12]));e=list(map(int,w[12:]))
    assert all(x%7==0 for x in e)
    A=[[sum((POW[(e[j]//7-3*H[r][j])%60] for j in range(11) if b[j]==k),K()) for k in range(5)] for r in range(1,5)]
    Am=[[sum(pow(2,(e[j]-21*H[r][j])%420,421) for j in range(11) if b[j]==k)%421 for k in range(5)] for r in range(1,5)]
    assert modular_rank(Am)==3
    N=[[K([Q(x) for x in v]) for v in row] for row in rec['kernel_basis']]
    free=rec['free_columns'];assert all(N[f][j]==int(i==j) for i,f in enumerate(free) for j in range(2))
    assert all(sum((A[r][j]*N[j][t] for j in range(5)),K())==0 for r in range(4) for t in range(2))
    if rec['exclusion']=='forced_zero_cycle_phase':assert not any(N[rec['column']])
    else:
        assert rec['exclusion']=='inconsistent_linear_norm_equations'
        eq=[[a.conj()*b,a*b.conj(),1-a*a.conj()-b*b.conj()] for a,b in N]
        witness=[K([Q(x) for x in v]) for v in rec['inconsistency_combination']]
        identity=[sum((witness[i]*eq[i][j] for i in range(5)),K()) for j in range(3)]
        assert not identity[0] and not identity[1] and identity[2]
report={'seed_exact_orthogonality':True,'ten_D_cycle_covers':True,'fourteen_certificates_independently_verified':True,'exclusions':dict(Counter(r['exclusion'] for r in records))}
Path('verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
