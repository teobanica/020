"""Independent exact verification of DDD supports, phases and root labels."""
from pathlib import Path
from itertools import product,combinations
from math import prod,gcd,lcm
from collections import Counter
import json,time
from supports import BASE,partitions,equations,row,offsets,orbits
from prepare_five import variables
from five_verify import forms,collision
ROOT=Path(__file__).resolve().parent

def mm(A,B):
    cols=list(zip(*B))
    return [[sum(x*y for x,y in zip(a,b)) for b in cols] for a in A]

def det(A):
    a=[r[:] for r in A];n=len(a);old=1;sign=1
    for k in range(n-1):
        if not a[k][k]:
            j=next((i for i in range(k+1,n) if a[i][k]),None)
            if j is None:return 0
            a[k],a[j]=a[j],a[k];sign=-sign
        z=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=z*a[i][j]-a[i][k]*a[k][j]
                assert num%old==0
                a[i][j]=num//old
        for i in range(k+1,n):a[i][k]=0
        old=z
    return sign*a[-1][-1]

def exact_family(Q,R,s,branch,labels):
    A,g=equations(Q,R);V=s['V'];U=s['U'];rank=s['rank'];diag=s['diagonal']
    ug=[sum(a*b for a,b in zip(r,g)) for r in U]
    theta=[sum(V[i][j]*(20//diag[j])*(ug[j]+2*branch[j]) for j in range(rank))%40 for i in range(6)]
    seed=[];basis=[]
    for P,off,base in ((BASE,0,0),(Q,3,8)):
        D=row(P,off);eps=offsets(P);lab=variables(P,base)
        seed.append([(8*labels[lab[j]]+20*eps[j]+sum(a*b for a,b in zip(D[j],theta)))%40 for j in range(11)])
        basis.append([[sum(D[j][i]*V[i][k] for i in range(6)) for k in range(rank,6)] for j in range(11)])
    ws=[(a-b)%40 for a,b in zip(*seed)];wb=[[a-b for a,b in zip(x,y)] for x,y in zip(*basis)]
    for const,coeff,P in zip(seed+[ws],basis+[wb],(BASE,Q,R)):
        for block in P:
            j=block[0];p=len(block)
            assert sorted((const[k]-const[j])%40 for k in block)==[40*k//p for k in range(p)]
            assert all(coeff[k]==coeff[j] for k in block)
    return {'modulus':40,'exponents':seed,'basis':basis,'supports':[BASE,Q,R]}

def extension_example():
    # Conference matrix of order six from the quadratic character modulo five.
    chi=lambda x:0 if x%5==0 else (1 if x%5 in (1,4) else -1)
    C=[[0 if i==j else (1 if i==0 or j==0 else chi(i-j)) for j in range(6)] for i in range(6)]
    assert C==list(map(list,zip(*C)))
    assert mm(C,C)==[[5*int(i==j) for j in range(6)] for i in range(6)]
    # H=I+iC; encode exponents of i, then dephase columns using the first row.
    raw=[[0 if i==j else (1 if C[i][j]==1 else 3) for j in range(6)] for i in range(6)]
    H=[[(x-y)%4 for x,y in zip(r,raw[0])] for r in raw]
    H=[[(x-r[0])%4 for x in r] for r in H]
    # Five rows: F_5 joined to the first five rows of H_6, in twentieth roots.
    E=[[(4*i*j)%20 for j in range(5)]+[(5*x)%20 for x in H[i]] for i in range(5)]
    certificates=[]
    for i,j in combinations(range(5),2):
        d=[(a-b)%20 for a,b in zip(E[i],E[j])]
        assert sorted(d[:5])==list(range(0,20,4))
        left=list(range(5,11));pairs=[]
        while left:
            a=left.pop(0);b=next(k for k in left if (d[k]-d[a])%20==10)
            left.remove(b);pairs.append([a,b])
        certificates.append({'rows':[i,j],'pentagon':list(range(5)),'pairs':pairs})
    return {'modulus':20,'exponents':E,'conference_matrix':C,'quaternary_H6_exponents':H,'quotient_certificates':certificates,
        'note':'Explicit regular 5x11 matrix with every quotient of type D. This is not a square completion.'}

def main():
    start=time.time();data=json.loads((ROOT/'supports.json').read_text());parts=partitions();reps=orbits(parts)
    assert len(parts)==6930 and len(set(parts))==6930 and len(reps)==26
    assert [(tuple(tuple(b) for b in o['representative']),o['size']) for o in data['support_orbits']]==reps
    phase=json.loads((ROOT/'phases.json').read_text());five=json.loads((ROOT/'five.json').read_text())
    assert len(phase)==180180
    five={(r['vi'],r['wi']):r for r in five};assert len(five)==15968
    stats=Counter();dimensions=Counter();supportdims=Counter();smithtypes=Counter();examples={};seedorders=Counter();no_binary_half=True
    for idx,s in enumerate(phase):
        vi,wi=divmod(idx,6930);assert (s['vi'],s['wi'])==(vi,wi)
        Q=data['support_orbits'][vi]['representative'];R=parts[wi]
        A,g=equations(Q,R);U=s['U'];V=s['V'];rank=s['rank'];diag=s['diagonal']
        assert abs(det(U))==abs(det(V))==1
        D=mm(mm(U,A),V)
        assert D==[[diag[i] if i==j and i<rank else 0 for j in range(6)] for i in range(7)]
        assert all(x>0 for x in diag) and all(b%a==0 for a,b in zip(diag,diag[1:]))
        assert all(x in (1,2,3,4) for x in diag)
        ug=[sum(a*b for a,b in zip(r,g)) for r in U]
        binary=all(x%2==0 for x in ug[rank:]);assert binary==s['binary_feasible']
        if not binary:
            stats['binary_excluded']+=1;assert (vi,wi) not in five;continue
        stats['binary_feasible']+=1
        assert all(x in (1,2,4) for x in diag) and rank>=3
        # A half-integral solution exists exactly when d_i*y_i=Ug_i mod 2
        # is consistent after the same unimodular coordinate changes.
        half_ok=all(d%2 or ug[i]%2==0 for i,d in enumerate(diag))
        no_binary_half &= not half_ok
        f=five[vi,wi];maps=[variables(P,b) for P,b in ((BASE,0),(Q,8),(R,16))];eq=list(zip(*maps))
        if not f['satisfiable']:
            assert 'collision' in f and collision(forms(eq))==f['collision']
            stats['five_linear_excluded']+=1;continue
        labels=f['labels'];assert len(labels)==24 and all(0<=x<5 for x in labels)
        assert all(sorted(labels[b:b+5])==list(range(5)) for b in (0,8,16))
        assert all((labels[a]-labels[b]-labels[c])%5==0 for a,b,c in eq)
        stats['feasible_supports']+=1;supportdims[6-rank]+=1;smithtypes[tuple(diag)]+=1
        for branch in product(*(range(d) for d in diag)):
            F=exact_family(Q,R,s,branch,labels);dimensions[6-rank]+=1
            order=lcm(*(40//gcd(40,(r[j]-r[0])%40) for r in F['exponents'] for j in range(11)))
            seedorders[order]+=1
            if order not in examples:
                div=40//order
                examples[order]={'vi':vi,'wi':wi,'branch':branch,'dimension':6-rank,
                    'modulus':order,'exponents':[[(x-r[0])%40//div for x in r] for r in F['exponents']],
                    'supports':F['supports'],'basis':F['basis']}
        if idx%10000==0:print('Verified',idx,'systems',flush=True)
    result={'verified':True,'partitions':len(parts),'relative_support_orbits':len(reps),'support_presentations':len(phase),**dict(stats),
        'feasible_supports_by_dimension':dict(sorted(supportdims.items())),
        'phase_components_by_dimension':dict(sorted(dimensions.items())),'phase_components':sum(dimensions.values()),
        'seed_orders_for_stored_witnesses':dict(sorted(seedorders.items())),
        'no_half_integral_binary_solutions':no_binary_half,
        'note':'Counts are redundant support/phase presentations, not matrix equivalence classes. All valid quintic labels are allowed; only one witness per support is stored.',
        'seconds':round(time.time()-start,3)}
    (ROOT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'examples.json').write_text(json.dumps(examples,indent=2)+'\n')
    (ROOT/'extension_5x11.json').write_text(json.dumps(extension_example(),indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
