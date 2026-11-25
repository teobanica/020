"""Independent checks of phase coverage, lifting certificates, and full families.

Uses no torus solver or finite-field search. All arithmetic is exact.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product, combinations, permutations
from collections import Counter, defaultdict
from math import gcd
import json
from collapse import BASE,partitions,row,equations
from prepare_lift import column_variables
from linear_check import normal_forms,collision

ROOT=Path(__file__).resolve().parent

def det(a):
    n=len(a)
    total=0
    for p in permutations(range(n)):
        sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        v=sign
        for i,j in enumerate(p):v*=a[i][j]
        total+=v
    return total

def primitive(B):
    d=len(B[0]);g=0
    if not d:return True
    for pick in combinations(range(4),d):g=gcd(g,det([B[j] for j in pick]))
    return g==1

def family(vi,wi,rec,labels,data,parts):
    P=data['support_orbits'][vi]['representative'];R=parts[wi]
    c=list(map(Q,rec['c']));B=rec['B'];d=rec['dimension']
    answer=[]
    for partition,offset,base in [(BASE,0,0),(P,2,9)]:
        c0,B0=row(partition,offset);var=column_variables(partition,base)
        const=[(c0[j]+3*sum(B0[j][i]*c[i] for i in range(4))+Q(labels[var[j]],7))%1 for j in range(11)]
        basis=[[sum(B0[j][i]*B[i][k] for i in range(4)) for k in range(d)] for j in range(11)]
        assert all((28*x).denominator==1 for x in const)
        answer.append((const,basis))
    u,v=answer
    w=([(a-b)%1 for a,b in zip(u[0],v[0])],[[a-b for a,b in zip(x,y)] for x,y in zip(u[1],v[1])])
    for (const,basis),partition in zip([u,v,w],[BASE,P,R]):
        for block in partition:
            a=block[0];p=len(block)
            assert sorted((const[j]-const[a])%1 for j in block)==[Q(k,p) for k in range(p)]
            assert all(basis[j]==basis[a] for j in block)
    return {'modulus':28,'exponents':[[int(28*x) for x in r[0]] for r in answer],
            'basis':[r[1] for r in answer],'supports':[BASE,P,R]}

def main():
    data=json.loads((ROOT/'collapsed.json').read_text());parts=partitions()
    grid=json.loads((ROOT/'grid_check.json').read_text())
    groups=defaultdict(list)
    for r in data['records']:groups[r['vi'],r['wi']].append(r)
    for check in grid:
        key=check['vi'],check['wi'];covered=set()
        A,b=equations(data['support_orbits'][key[0]]['representative'],parts[key[1]])
        for r in groups[key]:
            c=list(map(Q,r['c']));B=r['B'];dim=r['dimension']
            assert dim==4-check['rank'] and primitive(B)
            assert all((sum(x*y for x,y in zip(a,c))-rhs)%1==0 for a,rhs in zip(A,b))
            assert all(sum(a[i]*B[i][k] for i in range(4))==0 for a in A for k in range(dim))
            points={tuple(int(8*c[i])+sum(B[i][j]*t[j] for j in range(dim)) for i in range(4)) for t in product(range(8),repeat=dim)}
            points={tuple(x%8 for x in t) for t in points}
            assert len(points)==8**dim and not covered&points
            covered|=points
        assert len(covered)==check['grid_count'],(key,len(covered),check)
    lifts=json.loads((ROOT/'lifts.json').read_text());excluded=json.loads((ROOT/'excluded_lifts.json').read_text())
    assert len(lifts)==2571 and len(excluded)==424
    known=set();stats=Counter();sample=None;family_count=0
    for r in lifts:
        assert r['satisfiable']
        key=r['vi'],r['wi'];known.add(key);P=data['support_orbits'][key[0]]['representative'];R=parts[key[1]]
        labels=r['labels'];assert all(0<=x<7 for x in labels)
        assert all(sorted(labels[b:b+7])==list(range(7)) for b in [0,9,18])
        u,v,w=[column_variables(P,b) for P,b in [(BASE,0),(P,9),(R,18)]]
        assert all((labels[a]-labels[b]-labels[c])%7==0 for a,b,c in zip(u,v,w))
        for rec in groups[key]:
            F=family(*key,rec,labels,data,parts);family_count+=1;stats[rec['dimension']]+=1
            if sample is None and rec['dimension']==0 and any(Q(c).denominator==4 for c in rec['c']):
                sample={'vi':key[0],'wi':key[1],**F}
    for r in excluded:
        key=r['vi'],r['wi'];assert key not in known;known.add(key)
        P=data['support_orbits'][key[0]]['representative'];R=parts[key[1]]
        u,v,w=[column_variables(P,b) for P,b in [(BASE,0),(P,9),(R,18)]]
        col=min(set(BASE[0])&set(P[0])&set(R[0]));other=next(j for j in BASE[0] if j!=col)
        assert collision(normal_forms(list(zip(u,v,w)),[(u[col],0),(v[col],0),(u[other],1)]))==r['collision']
    assert known=={k for k,v in groups.items() if v}
    (ROOT/'example.json').write_text(json.dumps(sample,indent=2)+'\n')
    report={'verified':True,'grid_systems':len(grid),'max_smith_invariant':max(x['max_smith_invariant'] for x in grid),
            'liftable_support_pairs':len(lifts),'excluded_by_linear_collisions':len(excluded),
            'phase_components_with_lifts':family_count,'by_dimension':dict(sorted(stats.items())),
            'note':'Counts are redundant presentations, not monomial-equivalence classes. One root-label witness per support pair; general form allows every valid label assignment.'}
    (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    print('Example:',json.dumps(sample))

if __name__=='__main__':main()
