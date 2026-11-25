#!/usr/bin/env python3
"""Exact verification of the two maximal-dimensional FFF extension tests."""
from collections import Counter
from functools import lru_cache
from fractions import Fraction
from itertools import product
from math import factorial, gcd, lcm
from pathlib import Path
import json
import re
import sys
ROOT=Path(__file__).resolve().parent
M=210
PROFILES={'B':(7,2,2),'C':(5,3,3),'D':(5,2,2,2),'E':(3,3,3,2),'F':(3,2,2,2,2)}
PRESENTATIONS={'B':712800,'C':443520,'D':166320,'E':123200,'F':34650}

def records(name):
    return [tuple(map(int,s.split())) for s in (ROOT/name).read_text().splitlines()]

def monomial(a=0,b=0,s=0,t=0,p=0,q=0,e=0,group=0):
    assert 0<=group<=4
    return (a,b,s,t,p,q,*(int(group==i) for i in range(1,5)),e%M)

def quotient(x,y):
    return tuple(tuple(a-b for a,b in zip(u[:-1],v[:-1]))+((u[-1]-v[-1])%M,)
                 for u,v in zip(x,y))

def decode(key):
    assert len(key)==88
    row=tuple(monomial(a,b,s,t,p,q,e,g) for g,a,b,s,t,p,q,e in
              (key[i:i+8] for i in range(0,88,8)))
    assert row[0]==monomial()
    return row

def cover(row, profile=None):
    """Find literal prime cycles using monomial equality, independently of C++."""
    row = tuple(row)
    positions = {}
    for j, x in enumerate(row):
        positions.setdefault(x, []).append(j)
    wanted = tuple(sorted(profile or (), reverse=True))

    @lru_cache(None)
    def visit(mask, sizes):
        if not mask:
            return () if not sizes else None
        j = (mask & -mask).bit_length() - 1
        x = row[j]
        primes = sorted(set(sizes), reverse=True) if profile else (7, 5, 3, 2)
        for p in primes:
            if p > mask.bit_count():
                continue
            choices = [[k for k in positions.get(x[:-1] + ((x[-1] + M // p * a) % M,), ())
                        if mask >> k & 1] for a in range(1, p)]
            for others in product(*choices):
                block = (j,) + others
                assert len(set(block)) == p
                next_sizes = list(sizes)
                if profile:
                    next_sizes.remove(p)
                result = visit(mask ^ sum(1 << k for k in block), tuple(next_sizes))
                if result is not None:
                    return (block,) + result
        return None

    return visit((1 << len(row)) - 1, wanted)

def fixed_cover(row,blocks):
    assert sorted(j for b in blocks for j in b)==list(range(11))
    for b in blocks:
        assert len(b) in (2,3)
        first=row[b[0]]
        assert all(row[j][:-1]==first[:-1] for j in b)
        assert {(row[j][-1]-first[-1])%M for j in b}=={M//len(b)*k for k in range(len(b))}

def mm(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]

def determinant(a):
    a=[[Fraction(x) for x in row] for row in a];value=Fraction(1)
    for i in range(len(a)):
        pivot=next((j for j in range(i,len(a)) if a[j][i]),None)
        if pivot is None:return 0
        if pivot!=i:a[i],a[pivot]=a[pivot],a[i];value=-value
        value*=a[i][i]
        for j in range(i+1,len(a)):
            q=a[j][i]/a[i][i]
            a[j]=[x-q*y for x,y in zip(a[j],a[i])]
    return value

def check_original_components(data):
    # Reconstruct the phase equations independently in Python, then verify
    # the supplied unimodular identities and every orientation/component.
    parts=data['supports'];maps=[]
    for part in parts:
        maps.append({j:(i,k) for i,b in enumerate(part) for k,j in enumerate(b)})
    W=[];off=[];f=[]
    for j in range(11):
        (p,pp),(q,qp),(r,rp)=[m[j] for m in maps]
        row=[0]*8
        if p:row[p-1]=1
        if q:row[q+3]=-1
        W.append(row)
        off.append([-pp if p==0 else 0,qp if q==0 else 0,rp if r==0 else 0])
        f.append((rp if r else 0)-(pp if p else 0)+(qp if q else 0))
    A=[];B=[];rhs=[]
    for block in parts[2]:
        first=block[0]
        for j in block[1:]:
            A.append([x-y for x,y in zip(W[j],W[first])])
            B.append([x-y for x,y in zip(off[j],off[first])])
            rhs.append(f[j]-f[first])
    pc=data['phase_certificate'];U=pc['U'];V=pc['V']
    assert abs(determinant(U))==abs(determinant(V))==1
    assert mm(mm(U,A),V)==[[int(i==j and i<2) for j in range(8)] for i in range(6)]
    uf=[r[0] for r in mm(U,[[x] for x in rhs])];ub=mm(U,B)
    assert all(x%2==0 for x in uf[2:])
    allowed=[]
    for mask in range(8):
        sig=[1 if mask>>i&1 else -1 for i in range(3)]
        z=[sum(x*y for x,y in zip(row,sig)) for row in ub]
        if any(x%3 for x in z[2:]):continue
        allowed.append(mask)
        theta=[sum(row[k]*(3*uf[k]+2*z[k]) for k in range(2))%6 for row in V]
        rows=[]
        for h,start in ((0,0),(1,4)):
            row=[]
            for j in range(11):
                b,k=maps[h][j]
                powers=V[start+b-1][2:] if b else [0]*6
                e=3*k+theta[start+b-1] if b else 2*k*sig[h]
                row.append(monomial(*powers,e=35*e))
            rows.append(tuple(row))
        for row,blocks in zip((*rows,quotient(*rows)),parts):fixed_cover(row,blocks)
    assert allowed==pc['orientations'] and len(allowed)==2

def normal_form_check(u,thirds):
    base=[[0,1,2],[3,4],[5,6],[7,8],[9,10]]
    cross=[[0,1,2],[3,5],[4,6],[7,9],[8,10]]
    q1=[[0,1,2],[3,5],[4,6],[7,8],[9,10]]
    r1=[[0,1,2],[3,4],[5,6],[7,9],[8,10]]
    expected=[[base,cross,base],[base,q1,r1]]
    data=json.loads((ROOT/'normal_forms.json').read_text())
    assert len(data)==2
    assert [d['class_index'] for d in data]==[19818,19819]
    for h,d in enumerate(data):
        check_original_components(d)
        c=d['class_record'];pc=d['phase_certificate']
        assert c['dimension']==6 and c['orientations']==2 and c['torsion']==1
        assert pc['rank']==2 and pc['diagonal']==[1,1] and len(pc['orientations'])==2
        p,q,r=d['supports']
        if h==0: transformed=[p,r,q]
        else:
            perm=[0,1,2,7,8,9,10,3,4,5,6]
            transformed=[sorted([sorted(perm[j] for j in b) for b in part],
                                key=lambda b:(-len(b),b)) for part in (p,q,r)]
        assert transformed==expected[h]
        v=thirds[h]
        for row,blocks in zip((u,v,quotient(u,v)),expected[h]):fixed_cover(row,blocks)
        recovered=[quotient((u[5],),(u[3],))[0],quotient((u[9],),(u[7],))[0],u[3],v[3],u[7],v[7]]
        for i,value in enumerate(recovered):
            powers=[int(j==i) for j in range(6)]
            assert value==monomial(*powers)
    return {'support_types':2,'orientations_per_type':2,'phase_components':4,
            'parameters':6,'classification_source':'original fff census'}

def A_exceptional_equations(u,thirds):
    # In any A-containing triple every non-A quotient occupies at most two
    # mu_66 cosets. The three entries 1,s,p (or 1,t,q) must have a ratio in mu_66.
    assert [u[i] for i in (0,3,7)]==[monomial(),monomial(s=1),monomial(p=1)]
    for v in thirds:
        assert [v[i] for i in (0,3,7)]==[monomial(),monomial(t=1),monomial(q=1)]
    characters=[]
    for i,j in ((2,4),(3,5)):
        for first,second in ((i,None),(j,None),(i,j)):
            row=[0]*6;row[first]=1
            if second is not None:row[second]=-1
            characters.append(tuple(row))
    result={c+(35*d,) for c in characters for d in range(66)}
    assert len(result)==396
    return result

def main():
    if sys.flags.optimize:raise RuntimeError('Run without -O.')
    one=(monomial(),)*11
    u=tuple(monomial(e=70*j) for j in range(3))
    u+=(monomial(s=1),monomial(s=1,e=105),monomial(a=1,s=1),monomial(a=1,s=1,e=105))
    u+=(monomial(p=1),monomial(p=1,e=105),monomial(b=1,p=1),monomial(b=1,p=1,e=105))
    v0=tuple(monomial(e=140*j) for j in range(3))
    v0+=(monomial(t=1),monomial(t=1),monomial(t=1,e=105),monomial(t=1,e=105))
    v0+=(monomial(q=1),monomial(q=1),monomial(q=1,e=105),monomial(q=1,e=105))
    v1=v0[:7]+(monomial(q=1),monomial(q=1,e=105),monomial(b=1,q=1,e=105),monomial(b=1,q=1))
    thirds=(v0,v1)
    assert records('pair.txt')==[x[:6]+(x[-1],) for x in u]
    assert records('thirds.txt')==[x[:6]+(x[-1],) for v in thirds for x in v]
    normal=normal_form_check(u,thirds)
    a_ex=A_exceptional_equations(u,thirds)
    summary={};scan_ex=set();max_coefficient=0
    for p,sizes in PROFILES.items():
        count=factorial(11)
        for x in sizes:count//=x
        for n in Counter(sizes).values():count//=factorial(n)
        assert count==PRESENTATIONS[p]
        scans=[]
        for mode in ('primary','audit'):
            name=f'{mode}_{p}'
            pairs=records(name+'.txt.pairs');finals=records(name+'.txt')
            assert len(pairs)==len(set(pairs)) and len(finals)==len(set(finals))
            log=(ROOT/(name+'.log')).read_text()
            match=re.search(r'^tests (\d+) pairs (\d+) final (\d+) nodes ',log,re.M)
            assert match is not None
            assert tuple(map(int,match.groups()))==(count,len(pairs),len(finals))
            per_head=dict((int(a),int(b)) for a,b in re.findall(r'^family (\d+) fourth_rows (\d+)$',log,re.M))
            assert per_head=={h:sum(k[0]==h for k in finals) for h in range(2)}
            scans.append((set(pairs),set(finals)))
        assert scans[0]==scans[1]
        pairs,finals=scans[0]
        assert not finals
        assert len(pairs)==(61350 if p=='F' else 0)
        for key in pairs:
            row=decode(key)
            assert cover(row,sizes) is not None
            assert cover(quotient(row,u)) is not None
        for key in finals:
            h=key[0];assert h in (0,1);row=decode(key[1:])
            for v in (one,u,thirds[h]):assert cover(quotient(row,v)) is not None
        ex_counts=[]
        for stage in (0,1):
            ex=records(f'audit_{p}.txt.exceptional{stage}')
            assert len(ex)==len(set(ex))
            for row in ex:
                assert len(row)==7 and row[:6]!=(0,)*6
                assert next(x for x in row[:6] if x)>0 and 0<=row[-1]<M
                max_coefficient=max(max_coefficient,*map(abs,row[:6]))
                scan_ex.add(row[:6]+(11*row[-1],))
            ex_counts.append(len(ex))
        summary[p]={'presentations':count,'pair_records':len(pairs),
                    'triple_records':len(finals),'per_family':[sum(k[0]==h for k in finals) for h in range(2)],
                    'exceptional_counts':ex_counts}
        print(p,summary[p],flush=True)
    exception_set=scan_ex|a_ex
    assert len(scan_ex)==21043 and len(exception_set)==21403 and max_coefficient==4
    order=83;parameter_exponents=(1,7,51,15,73,20);level=6*order
    for row in exception_set:
        dot=sum(a*b for a,b in zip(row[:6],parameter_exponents))
        assert (2310*dot-order*row[-1])%(2310*order)!=0
    seeds=[]
    for row in (u,*thirds):
        assert all(x[-1]%35==0 for x in row)
        seeds.append(tuple((order*(x[-1]//35)+6*sum(a*b for a,b in
                       zip(x[:6],parameter_exponents)))%level for x in row))
    for h in range(2):
        eu,ev=seeds[0],seeds[h+1]
        for symbolic,actual in ((u,eu),(thirds[h],ev),
                               (quotient(u,thirds[h]),tuple((a-b)%level for a,b in zip(eu,ev)))):
            blocks=cover(symbolic,PROFILES['F']);assert blocks is not None
            for block in blocks:
                assert {(actual[j]-actual[block[0]])%level for j in block}=={
                    level//len(block)*k for k in range(len(block))}
        assert eu[0]==ev[0]==0
        assert lcm(*(level//gcd(e,level) for row in (eu,ev) for e in row))==level
    examples={'parameter_order':order,'parameter_exponents':parameter_exponents,
              'minimal_Butson_level':level,'u':seeds[0],'v0':seeds[1],'v1':seeds[2],
              'exceptional_equations_satisfied':0}
    summary={'profiles':summary,'normal_form':normal,'A_exceptional_equations':len(a_ex),
             'scan_exceptional_equations':len(scan_ex),'exceptional_union':len(exception_set),
             'root_denominator':2310,'max_exceptional_coefficient':max_coefficient,
             'maximal_three_row_examples':examples,'passed':True}
    (ROOT/'exceptional_union.txt').write_text(''.join(' '.join(map(str,r))+'\n' for r in sorted(exception_set)))
    (ROOT/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
