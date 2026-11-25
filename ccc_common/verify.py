#!/usr/bin/env python3
"""Independent exact verification of the common-block CCC breadth test.

The two complete searches establish exclusion; this checker validates their
inputs, census accounting, agreement, positive row certificates and the finite
exceptional sets. It uses literal Laurent monomials, not phase graphs.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations, product
from math import factorial
from pathlib import Path
import json
import re
import sys

ROOT=Path(__file__).resolve().parent
M=210
PROFILES={'B':(7,2,2),'C':(5,3,3),'D':(5,2,2,2),
          'E':(3,3,3,2),'F':(3,2,2,2,2)}
PRESENTATIONS={'B':712800,'C':443520,'D':166320,'E':123200,'F':34650}
COUNTS={'B':(0,0),'C':(4383,852),'D':(672,426),'E':(0,0),'F':(0,0)}


def records(name):
    return [tuple(map(int,s.split())) for s in (ROOT/name).read_text().splitlines()]


def monomial(s=0,r=0,t=0,b=0,e=0,group=0):
    assert 0<=group<=4
    return (s,r,t,b,*(int(group==i) for i in range(1,5)),e%M)


def quotient(x,y):
    return tuple(tuple(a-b for a,b in zip(u[:-1],v[:-1]))+((u[-1]-v[-1])%M,)
                 for u,v in zip(x,y))


def decode(key):
    assert len(key)==66
    row=tuple(monomial(s,r,t,b,e,g) for g,s,r,t,b,e in
              (key[i:i+6] for i in range(0,66,6)))
    assert row[0]==monomial()
    return row


def circle(k,h,group=1):
    row=tuple(monomial(e=42*k*j) for j in range(5))
    row+=tuple(monomial(e=70*h*j,group=group) for j in range(3))
    powers=((0,0,0,0),(-1,1,0,0),(0,0,-1,1))[h]
    row+=tuple(monomial(*powers,e=105+70*h*j,group=group) for j in range(3))
    return row


def divide(p,q):
    p=list(p);answer=[0]*(len(p)-len(q)+1)
    while len(p)>=len(q):
        c=p[-1];i=len(p)-len(q);answer[i]=c
        for j,x in enumerate(q):p[i+j]-=c*x
        while p and p[-1]==0:p.pop()
    assert not p
    return tuple(answer)


@lru_cache(None)
def cyclotomic(n):
    p=(-1,)+(0,)*(n-1)+(1,)
    for d in range(1,n):
        if n%d==0:p=divide(p,cyclotomic(d))
    return p


@lru_cache(None)
def root_coefficient(e):
    p=[0]*M;p[e%M]=1;q=cyclotomic(M);d=len(q)-1
    for k in range(M-1,d-1,-1):
        c=p[k]
        for j,x in enumerate(q):p[k-d+j]-=c*x
    return tuple(p[:d])


def inner_product(x,y):
    result={}
    for v in quotient(x,y):
        value=root_coefficient(v[-1])
        result[v[:-1]]=tuple(a+b for a,b in zip(result.get(v[:-1],(0,)*len(value)),value))
    return {k:v for k,v in result.items() if any(v)}


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



def main():
    if sys.flags.optimize:
        raise RuntimeError('Run without -O: assertions are part of this checker.')
    heads=[(0,)+f for f in permutations(range(1,5))
           if len({(v-j)%5 for j,v in enumerate((0,)+f)})==5]
    assert heads==records('heads.txt')
    assert heads==[tuple(a*j%5 for j in range(5)) for a in (2,3,4)]
    # Rank zero means every R block lies in one P/Q intersection cell.
    # Among all 4620 Q partitions, only Q=P has cells fitting 5,3,3.
    base=(frozenset(range(5)),frozenset(range(5,8)),frozenset(range(8,11)))
    zero_rank=[];partition_count=0
    for big in combinations(range(11),5):
        rest=sorted(set(range(11))-set(big))
        for two in combinations(rest[1:],2):
            tri=frozenset((rest[0],)+two)
            q=(frozenset(big),tri,frozenset(set(rest)-tri))
            partition_count+=1
            cells=[p&v for p in base for v in q if p&v]
            if sorted(map(len,cells))==[3,3,5]:zero_rank.append(q)
    assert partition_count==4620 and zero_rank==[base]
    one=(monomial(),)*11
    u=tuple(monomial(e=42*j) for j in range(5))
    u+=tuple(monomial(s=1,e=70*j) for j in range(3))
    u+=tuple(monomial(r=1,e=70*j) for j in range(3))
    assert records('pair.txt')==[x[:4]+(x[-1],) for x in u]
    old=[]
    for head in heads:
        v=tuple(monomial(e=42*j) for j in head)
        v+=tuple(monomial(t=1,e=140*j) for j in range(3))
        v+=tuple(monomial(b=1,e=140*j) for j in range(3))
        old.append(v)
        for x,y in ((u,one),(v,one),(u,v)):
            assert cover(quotient(x,y),PROFILES['C']) is not None
    summary={};exception_set=set();max_coefficient=0;circles=set();fixed=Counter()
    for p,sizes in PROFILES.items():
        count=factorial(11)
        for x in sizes:count//=x
        for m in Counter(sizes).values():count//=factorial(m)
        assert count==PRESENTATIONS[p]
        scans=[]
        for mode in ('primary','audit'):
            name=f'{mode}_{p}'
            pairs=records(name+'.txt.pairs');finals=records(name+'.txt')
            assert len(pairs)==len(set(pairs)) and len(finals)==len(set(finals))
            match=re.search(r'^tests (\d+) pairs (\d+) final (\d+) nodes ',
                            (ROOT/(name+'.log')).read_text(),re.M)
            assert match is not None
            assert tuple(map(int,match.groups()))==(count,len(pairs),len(finals))
            scans.append((set(pairs),set(finals)))
        assert scans[0]==scans[1]
        pairs,finals=scans[0]
        assert (len(pairs),len(finals))==COUNTS[p]
        if p in ('C','D'):
            assert Counter(row[0] for row in finals)==Counter({i:(284 if p=='C' else 142) for i in range(3)})
        for key in pairs:
            x=decode(key)
            assert cover(x,sizes) is not None
            assert cover(quotient(x,u)) is not None
        for rec in finals:
            i,key=rec[0],rec[1:]
            assert 0<=i<len(heads)
            x=decode(key)
            for y in (one,u,old[i]):
                assert cover(quotient(x,y)) is not None
            k=x[1][-1]//42
            a=(2,3,4)[i]
            assert k in set(range(5))-{0,1,a}
            assert x[:5]==tuple(monomial(e=42*k*j) for j in range(5))
            h=quotient(x[6:7],x[5:6])[0][-1]//70
            assert h in (0,1,2)
            template=circle(k,h)
            assert quotient(x[5:],(x[5],)*6)==quotient(template[5:],(template[5],)*6)
            assert p==('D' if h==0 else 'C')
            groups=tuple(key[j] for j in range(0,66,6))
            if max(groups):
                assert groups==(0,)*5+(1,)*6 and x[5]==monomial(group=1)
                circles.add((i,k,h))
            else:fixed[p]+=1
        ex_counts=[]
        for stage in (0,1):
            ex=records(f'audit_{p}.txt.exceptional{stage}')
            assert len(ex)==len(set(ex))
            for row in ex:
                assert len(row)==5 and row[:4]!=(0,)*4
                assert next(x for x in row[:4] if x)>0 and 0<=row[-1]<M
                max_coefficient=max(max_coefficient,*map(abs,row[:4]))
            exception_set.update(ex);ex_counts.append(len(ex))
        summary[p]={'presentations':count,'pair_records':len(pairs),
                    'triple_records':len(finals),'exceptional_counts':ex_counts}
        print(p,summary[p],flush=True)
    expected={(i,k,h) for i,a in enumerate((2,3,4))
              for k in set(range(5))-{0,1,a} for h in (0,1,2)}
    assert circles==expected
    assert fixed==Counter({'C':840,'D':420})
    assert len(exception_set)==3608 and max_coefficient==3
    compatibility={}
    zero=(0,)*8
    phase=(0,)*4+(1,-1,0,0)
    for i,a in enumerate((2,3,4)):
        vertices=sorted((k,h) for ii,k,h in circles if ii==i)
        edges=[]
        for k,h in vertices:
            x=circle(k,h)
            for y in (one,u,old[i]):
                assert cover(quotient(x,y)) is not None
                assert inner_product(x,y)=={}
        for (k,h),(l,j) in combinations_with_replacement(vertices,2):
            actual=inner_product(circle(k,h,1),circle(l,j,2))
            predicted={}
            if k==l:predicted[zero]=(5,)+(0,)*(len(cyclotomic(M))-2)
            if h==j:predicted[phase]=(6,)+(0,)*(len(cyclotomic(M))-2)
            assert actual==predicted
            if not actual:edges.append(((k,h),(l,j)))
        assert len(vertices)==6 and len(edges)==6
        degrees=Counter(v for edge in edges for v in edge)
        assert set(degrees.values())=={2}
        adjacency={frozenset(edge) for edge in edges}
        assert not any(all(frozenset(pair) in adjacency for pair in combinations(triple,2))
                       for triple in combinations(vertices,3))
        compatibility[a]={'vertices':vertices,'edges':edges,'maximum_clique':2}
    summary={'profiles':summary,'heads':heads,'support_partitions_checked':partition_count,
             'zero_rank_supports':len(zero_rank),'exceptional_union':len(exception_set),
             'max_exceptional_coefficient':max_coefficient,'circles':sorted(circles),
             'fixed_specializations':dict(fixed),'compatibility':compatibility,'passed':True}
    (ROOT/'exceptional_union.txt').write_text(''.join(' '.join(map(str,v))+'\n'
                                                   for v in sorted(exception_set)))
    (ROOT/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
