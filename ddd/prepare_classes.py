from pathlib import Path
from collections import Counter
from math import prod
import json
from supports import BASE,partitions
ROOT=Path(__file__).resolve().parent
def block_map(P):
    out=[None]*11
    for k,b in enumerate(P):
        for j in b:out[j]=k
    return out
def main():
    data=json.loads((ROOT/'supports.json').read_text());parts=partitions();P=block_map(BASE)
    maps=[block_map(r) for r in parts]
    qs=[block_map(o['representative']) for o in data['support_orbits']]
    phase=json.loads((ROOT/'phases.json').read_text())
    five={(r['vi'],r['wi']):r for r in json.loads((ROOT/'five.json').read_text())}
    tensors=set()
    with (ROOT/'classes_input.txt').open('w') as f:
        f.write(str(len(phase))+'\n')
        for s in phase:
            vi,wi=s['vi'],s['wi'];Q,R=qs[vi],maps[wi]
            cols=[16*a+4*b+c for a,b,c in zip(P,Q,R)]
            tensors.add(tuple(sorted(cols)))
            ok=s['binary_feasible'] and five[vi,wi]['satisfiable']
            f.write(' '.join(map(str,[vi,wi,6-s['rank'],prod(s['diagonal']),int(ok),*cols]))+'\n')
    print('Distinct raw intersection tensors:',len(tensors),flush=True)
if __name__=='__main__':main()
