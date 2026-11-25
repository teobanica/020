#!/usr/bin/env python3
"""Independent exact verifier for the four-parameter EEE extension test."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial, gcd, lcm
from pathlib import Path
import json
import re
import sys

ROOT=Path(__file__).resolve().parent
M=420
PROFILES={'B':(7,2,2),'C':(5,3,3),'D':(5,2,2,2),
          'E':(3,3,3,2),'F':(3,2,2,2,2)}
PRESENTATIONS={'B':712800,'C':443520,'D':166320,'E':123200,'F':34650}


def records(name):
    return [tuple(map(int,s.split())) for s in (ROOT/name).read_text().splitlines()]


def monomial(s=0,t=0,a=0,b=0,e=0,group=0):
    assert 0<=group<=4
    return (s,t,a,b,*(int(group==i) for i in range(1,5)),e%M)


def quotient(x,y):
    return tuple(tuple(a-b for a,b in zip(u[:-1],v[:-1]))+((u[-1]-v[-1])%M,)
                 for u,v in zip(x,y))


def decode(key):
    assert len(key)==66
    row=tuple(monomial(s,t,a,b,e,g) for g,s,t,a,b,e in
              (key[i:i+6] for i in range(0,66,6)))
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




def check_A_graphs():
    # All labelled simple bipartite graphs for AE(E/F), with E on the left.
    result={}
    for name,right in [('EE',(3,3,3,2)),('EF',(3,2,2,2,2))]:
        count=0;n=len(right)
        choices=[list(combinations(range(n),d)) for d in (3,3,3,2)]
        for rows in product(*choices):
            if tuple(sum(j in row for row in rows) for j in range(n))!=right:continue
            count+=1;adj={i:set() for i in range(4+n)}
            for i,row in enumerate(rows):
                for j in row:adj[i].add(4+j);adj[4+j].add(i)
            seen={0};todo=[0]
            while todo:
                for j in adj[todo.pop()]-seen:seen.add(j);todo.append(j)
            assert len(seen)==4+n
        assert count>0
        result[name]=count
    return result


def normal_form_check():
    # For the unique dimension-four support type, the five-column core has
    # pair supports 01,02,12 and triangle supports 234,134,034.
    cores=[]
    for sigma,tau in product((1,-1),repeat=2):
        for x,y in product(range(12),repeat=2):
            u=(0,6,x,(x+4*sigma)%12,(x-4*sigma)%12)
            v=(0,y,6,(y+4*tau)%12,(y-4*tau)%12)
            q=tuple((a-b)%12 for a,b in zip(u,v))
            if (q[1]-q[2])%12!=6:continue
            if set(q[j] for j in (0,3,4))!={0,4,8}:continue
            cores.append((sigma,tau,x,y))
    expected=[(sigma,-sigma,x,x) for sigma in (1,-1) for x in (3,9)]
    assert sorted(cores)==sorted(expected)
    # Exchanging the last two core columns handles sigma=-1. Conjugation
    # exchanges x=3 and x=9, with that same column swap.
    B=((0,0,0,0,0),(0,6,3,7,11),(0,3,6,11,7))
    for sigma,tau,x,y in cores:
        row1=(0,6,x,(x+4*sigma)%12,(x-4*sigma)%12)
        row2=(0,y,6,(y+4*tau)%12,(y-4*tau)%12)
        found=False
        for conjugate,swap in product((False,True),repeat=2):
            cols=(0,1,2,4,3) if swap else (0,1,2,3,4)
            sign=-1 if conjugate else 1
            test=tuple(tuple(sign*row[j]%12 for j in cols) for row in B[1:])
            if test==(row1,row2):found=True
        assert found
    orientations=[]
    for sigma,tau in product((1,-1),repeat=2):
        q=tuple(monomial(e=140*(sigma-tau)*j) for j in range(3))
        if cover(q,(3,)) is not None:orientations.append((sigma,tau))
    assert orientations==[(1,-1),(-1,1)]
    return {'core_solutions':len(cores),
            'triangle_orientation_choices':len(orientations)**2,
            'total_components':len(cores)*len(orientations)**2}


def main():
    if sys.flags.optimize:raise RuntimeError('Run without -O.')
    one=(monomial(),)*11
    u=tuple(monomial(e=35*j) for j in (0,6,3,7,11))
    u+=tuple(monomial(s=1,e=140*j) for j in range(3))
    u+=tuple(monomial(a=1,e=140*j) for j in range(3))
    v=tuple(monomial(e=35*j) for j in (0,3,6,11,7))
    v+=tuple(monomial(t=1,e=280*j) for j in range(3))
    v+=tuple(monomial(b=1,e=280*j) for j in range(3))
    assert records('triple.txt')==[x[:4]+(x[-1],) for x in u+v]
    example=json.loads((ROOT/'eee_example4.json').read_text())
    assert example['dimension']==4 and example['modulus']==12
    assert example['supports']==[
        [[0,1],[2,3,4],[5,6,7],[8,9,10]],
        [[0,2],[1,3,4],[5,6,7],[8,9,10]],
        [[1,2],[0,3,4],[5,6,7],[8,9,10]]]
    for row,basis,expected in zip(example['exponents'],example['basis'],(u,v)):
        actual=tuple(monomial(*(powers[i] for i in (1,2,0,3)),e=35*e)
                     for e,powers in zip(row,basis))
        assert actual==expected
    for x,y in ((u,one),(v,one),(u,v)):
        assert cover(quotient(x,y),PROFILES['E']) is not None
    assert u[2][-1]==v[1][-1]==M//4 and 66%4!=0
    graphs=check_A_graphs();normal=normal_form_check()
    summary={};exception_set=set();max_coefficient=0
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
        assert not finals
        assert len(pairs)=={'B':0,'C':0,'D':0,'E':4184,'F':6418}[p]
        for key in pairs:
            x=decode(key)
            assert cover(x,sizes) is not None
            assert cover(quotient(x,u)) is not None
        for key in finals:
            x=decode(key)
            for y in (one,u,v):assert cover(quotient(x,y)) is not None
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
    assert len(exception_set)==5108 and max_coefficient==3
    # Explicit maximal three-row example: the old parameters are the
    # powers 1,2,4,7 of a primitive seventeenth root.
    parameter_exponents=(1,2,4,7)
    for r in exception_set:
        dot=sum(a*b for a,b in zip(r[:4],parameter_exponents))
        assert (420*dot-17*r[-1])%(420*17)!=0
    seed=[]
    for row in (u,v):
        assert all(x[-1]%35==0 for x in row)
        seed.append(tuple((17*(x[-1]//35)+
                           12*sum(a*b for a,b in zip(x[:4],parameter_exponents)))%204
                          for x in row))
    for row,blocks in zip((seed[0],seed[1],
                          tuple((a-b)%204 for a,b in zip(*seed))),example['supports']):
        for block in blocks:
            assert {(row[j]-row[block[0]])%204 for j in block}=={
                204//len(block)*k for k in range(len(block))}
    assert all(row[0]==0 for row in seed)
    assert lcm(*(204//gcd(e,204) for row in seed for e in row))==204
    witness={'parameter_order':17,'parameter_exponents':parameter_exponents,
             'minimal_Butson_level':204,'exponent_rows':seed,
             'exceptional_equations_satisfied':0}
    summary={'profiles':summary,'A_graph_counts':graphs,'normal_form':normal,
             'exceptional_union':len(exception_set),
             'max_exceptional_coefficient':max_coefficient,
             'maximal_three_row_example':witness,'passed':True}
    (ROOT/'exceptional_union.txt').write_text(''.join(' '.join(map(str,v))+'\n'
                                                   for v in sorted(exception_set)))
    (ROOT/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
