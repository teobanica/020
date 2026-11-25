#!/usr/bin/env python3
"""Check the two complete one-parameter DDD classes with fixed pair 21."""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path
import hashlib
import json
import sys

from normalize import canonical, U, CASE_COUNTS
from circle_bounds import check_sharp

ROOT = Path(__file__).resolve().parent
PAIR_HASH = '0d251854aedf2de521bbbbd3c587d6da247207c9db2e084a12c77ebc60f4a496'
L = 5040


def multiply(a, b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def determinant(a):
    a = [[Fraction(x) for x in row] for row in a]
    value = Fraction(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i,len(a)) if a[j][i]), None)
        if pivot is None:
            return 0
        if pivot != i:
            a[i],a[pivot] = a[pivot],a[i]
            value = -value
        value *= a[i][i]
        for j in range(i+1,len(a)):
            t = a[j][i]/a[i][i]
            a[j] = [x-t*y for x,y in zip(a[j],a[i])]
    return value


def check_labels(directory, class_id):
    c = json.loads((directory/f'class{class_id}.json').read_text())
    assert c['class_id']==class_id and c['dimension']==1
    supports = c['supports']
    positions = [{j:(i,k) for i,block in enumerate(part) for k,j in enumerate(block)}
                 for part in supports]
    weights = []
    signs = []
    for j in range(11):
        (p,pi),(q,qi),(r,ri) = [m[j] for m in positions]
        w = [0]*6
        if p: w[p-1] = 1
        if q: w[q+2] = -1
        weights.append(w)
        signs.append((ri if r else 0)-(pi if p else 0)+(qi if q else 0))
    edges = [(j,b[0]) for b in supports[2] for j in b[1:]]
    A = [[x-y for x,y in zip(weights[j],weights[h])] for j,h in edges]
    g = [signs[j]-signs[h] for j,h in edges]
    S,V = c['U'],c['V']
    diagonal = [1,1,1,1,2]
    assert A==c['A'] and g==c['g'] and c['diagonal']==diagonal
    assert abs(determinant(S))==abs(determinant(V))==1
    assert multiply(multiply(S,A),V)==[
        [diagonal[i] if i==j and i<5 else 0 for j in range(6)] for i in range(7)]
    sg = [sum(a*b for a,b in zip(row,g)) for row in S]
    assert all(x%2==0 for x in sg[5:])
    labels = [(0,)+p for p in permutations(range(1,5))]
    results = set()
    for pa,pb,pc in product(range(24),repeat=3):
        pent = []
        for part,k in zip(supports,(pa,pb,pc)):
            row = [0]*11
            for j,v in zip(part[0],labels[k]):row[j]=v
            pent.append(row)
        raw = [z-x+y for x,y,z in zip(*pent)]
        rhs = [raw[j]-raw[h] for j,h in edges]
        sf = [sum(x*y for x,y in zip(row,rhs)) for row in S]
        if any(x%5 for x in sf[5:]):
            continue
        for branch in range(2):
            y = [8*(pow(d,-1,5)*sf[i]%5)+(20//d)*(sg[i]+(2*branch if d==2 else 0))
                 for i,d in enumerate(diagonal)]
            phases = [sum(x*t for x,t in zip(row,y)) for row in V]
            exp = []
            slope = []
            for h in range(2):
                for j in range(11):
                    b,k = positions[h][j]
                    coord = 3*h+b-1
                    exp.append((8*pent[h][j]+(20*k+phases[coord] if b else 0))%40)
                    slope.append(V[coord][5] if b else 0)
            results.add((class_id,1,pa,pb,pc,branch,*exp,*slope))
    saved = [tuple(map(int,s.split())) for s in (directory/'labels.txt').read_text().splitlines()]
    assert len(results)==len(saved)==960 and results==set(saved)
    return results


def cover(row):
    # Entries are Laurent monomials in the old and new free phases, with
    # their last coordinate a root exponent modulo L.
    @lru_cache(None)
    def solve(mask):
        if not mask:
            return True
        h = (mask&-mask).bit_length()-1
        a,g,e = row[h]
        for p in (2,3,5,7):
            if p>mask.bit_count():continue
            choices = [[j for j in range(11) if mask>>j&1 and
                        row[j]==(a,g,(e+L//p*k)%L)] for k in range(1,p)]
            for other in product(*choices):
                used = (1<<h)+sum(1<<j for j in other)
                if solve(mask^used):return True
        return False
    return solve(2047)


def quotient(x,y):
    return tuple((a-b,g-h,(e-f)%L) for (a,g,e),(b,h,f) in zip(x,y))


def check(directory=ROOT, class_id=3203):
    assert not sys.flags.optimize
    for prime in (2521,3361):
        assert all(prime%d for d in range(2,int(prime**0.5)+1))
    labels = check_labels(directory,class_id)
    families = CASE_COUNTS[class_id]
    targets = [list(map(int,l.split())) for l in (directory/'target.txt').read_text().splitlines()]
    keys = [tuple(v for d,e in zip(t[12:],t[1:12]) for v in (d,e)) for t in targets]
    assert sorted(set(canonical(row) for row in labels))==keys and len(keys)==families
    assert [r[0] for r in targets]==list(range(families))
    equal = ((0,0,0),)*11
    u = tuple((0,0,126*e) for e in U)
    old = [tuple((d,0,e*126) for d,e in zip(row[12:],row[1:12])) for row in targets]
    for v in old:
        active=[j for j,(a,g,e) in enumerate(v) if a]
        assert len(active)==2 and U[active[0]]==U[active[1]]
        assert (v[active[0]][2]-v[active[1]][2])%L==L//2
        assert all(e%(2*126)==0 for a,g,e in v)
        assert all(cover(row) for row in (u,v,quotient(u,v)))
    pair = ROOT/'g21_pair_families.txt'
    assert hashlib.sha256(pair.read_bytes()).hexdigest()==PAIR_HASH
    counts = Counter()
    with pair.open() as f:
        for l in f:
            p,idx,n,*r = l.split()
            assert len(r)==22 and int(n)==1+max(map(int,r[:11]))<=3
            counts[p]+=1
    assert counts=={'C':1196,'D':51044,'E':17632,'F':269136}
    main = directory/'extensions.txt'
    audit = directory/'audit_extensions.txt'
    assert main.read_bytes()==audit.read_bytes()
    seen=set();points=defaultdict(set);types=Counter();group_sizes=set()
    for line in main.read_text().splitlines():
        w=tuple(map(int,line.split()))
        assert len(w)==35 and w not in seen
        seen.add(w);h,z=w[:2];g=w[2:13];a=w[13:24];e=w[24:]
        assert 0<=h<families and max(g)<=1 and min(g)==0
        if z>=0:
            assert 0<=z<L and z%84==0 and not any(a)
            v=tuple((0,0,(f+d*z)%L) for d,k,f in old[h]);points[h].add(z)
        else:
            assert z==-1
            v=old[h]
        x=tuple(zip(a,g,e))
        assert all(cover(quotient(x,y)) for y in (equal,u,v))
        types[(z==-1,max(g)+1)]+=1
        if max(g)==1:group_sizes.add(tuple(sorted(Counter(g).values())))
    expected_types = {3203:(4016,96,48256,384),3217:(2904,72,32400,288)}
    assert types==dict(zip([(True,1),(True,2),(False,1),(False,2)], expected_types[class_id]))
    domains={(h,q) for h in range(families) for q in {-1}|points[h]}
    bounds={}
    for name in ('graphs.txt','audit_graphs.txt'):
        rows=[tuple(map(int,s.split())) for s in (directory/name).read_text().splitlines()]
        assert len(rows)==len(domains)=={3203:1408,3217:512}[class_id] and {(r[0],r[1]) for r in rows}==domains
        assert all(0<=r[-1]<=3 for r in rows)
        bounds[name]=dict(sorted(Counter(r[-1] for r in rows).items()))
    summary={'support_class':class_id,'labelled_components':len(labels),'normalized_families':families,
             'pair_records':sum(counts.values()),'records_free_in_z':types[(True,1)]+types[(True,2)],
             'records_at_fixed_z':types[(False,1)]+types[(False,2)],'additional_phase_group_sizes':sorted(group_sizes),
             'exceptional_values_union':len(set.union(*points.values())),'graphs':len(domains),
             'clique_histograms':bounds,'uniform_total_row_bound':6,'passed':True}
    if class_id==3217:
        sharp=check_sharp(directory,old,u,cover,quotient)
        summary['uniform_total_row_bound']=sharp['maximum_total_rows']
        summary['sharp_profile']=sharp
    (directory/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(f'DDD class {class_id}: proved upper bound {summary["uniform_total_row_bound"]}, {len(domains)} graphs verified.')


if __name__=='__main__':
    check()
    check(ROOT/'3217',3217)
