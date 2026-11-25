from pathlib import Path
import json
from supports import BASE,partitions
ROOT=Path(__file__).resolve().parent

def variables(P,prime,base):
    out=[None]*11;k=base
    for block in P:
        if len(block)==prime:
            for j in block:out[j]=k;k+=1
        else:
            for j in block:out[j]=k
            k+=1
    assert k==base+7
    return out

def main():
    data=json.loads((ROOT/'supports.json').read_text());parts=partitions()
    u=variables(BASE,5,0)
    with (ROOT/'five_input.txt').open('w') as f:
        f.write(str(len(data['support_orbits'])*len(parts))+'\n')
        for vi,o in enumerate(data['support_orbits']):
            v=variables(o['representative'],5,7)
            for wi,R in enumerate(parts):
                w=variables(R,5,14)
                f.write(f'{vi} {wi}\n')
                for a,b,c in zip(u,v,w):f.write(f'{a} {b} {c}\n')
    print('Prepared quintic label systems.',flush=True)
if __name__=='__main__':main()
