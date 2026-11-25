"""Classify chosen CCC support triples under column and row permutations.

The 3x3x3 intersection tensor determines an ordered triple of set partitions
up to column permutation. Swap the two 3-blocks in each partition and permute
the three quotient partitions. This is not matrix-equivalence classification.
"""
from itertools import permutations,product
from functools import lru_cache
from pathlib import Path
from collections import Counter
from math import prod
import json
from supports import BASE,partitions
ROOT=Path(__file__).resolve().parent

def block_map(P):
    a=[None]*11
    for k,b in enumerate(P):
        for j in b:a[j]=k
    return a

def tensor(P,Q,R):
    a=[0]*27
    for x,y,z in zip(block_map(P),block_map(Q),block_map(R)):a[9*x+3*y+z]+=1
    return tuple(a)

TRANSFORMS=[]
for axis in permutations(range(3)):
    for bits in product(range(2),repeat=3):
        m=[]
        for p in product(range(3),repeat=3):
            q=[p[a] for a in axis]
            q=[3-x if b and x else x for x,b in zip(q,bits)]
            m.append(9*q[0]+3*q[1]+q[2])
        assert sorted(m)==list(range(27));TRANSFORMS.append(m)

@lru_cache(None)
def canonical(t):
    return min(tuple(t[j] for j in m) for m in TRANSFORMS)

def main():
    data=json.loads((ROOT/'supports.json').read_text());parts=partitions()
    five=json.loads((ROOT/'five.json').read_text());three=json.loads((ROOT/'three.json').read_text())
    classes={}
    for f,s in zip(five,three):
        Q=data['support_orbits'][f['vi']]['representative'];R=parts[f['wi']]
        c=canonical(tensor(BASE,Q,R));ok=f['satisfiable'] and bool(s['orientations'])
        if c in classes:
            assert classes[c]['feasible']==ok
            assert classes[c]['dimension']==4-s['rank']
            assert classes[c]['torsion_order']==prod(s['diagonal'])
        else:classes[c]={'tensor':c,'feasible':ok,'dimension':4-s['rank'],'torsion_order':prod(s['diagonal']),'representative':[f['vi'],f['wi']]}
    result={'total_support_classes':len(classes),'feasible_support_classes':sum(x['feasible'] for x in classes.values()),
        'feasible_by_dimension':dict(sorted(Counter(x['dimension'] for x in classes.values() if x['feasible']).items())),
        'feasible_by_torsion_order':dict(sorted(Counter(x['torsion_order'] for x in classes.values() if x['feasible']).items())),
        'classes':list(classes.values()),'scope':'Chosen support triples modulo S_11 column permutations, S_3 quotient permutations, and swapping equal-size blocks. Not matrix equivalence classes.'}
    (ROOT/'support_classes.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='classes'},indent=2),flush=True)

if __name__=='__main__':main()
