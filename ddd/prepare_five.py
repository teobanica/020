from pathlib import Path
import json
from supports import BASE,partitions
ROOT=Path(__file__).resolve().parent

def variables(P,base):
    out=[None]*11;k=base
    for block in P:
        if len(block)==5:
            for j in block:out[j]=k;k+=1
        else:
            for j in block:out[j]=k
            k+=1
    assert k==base+8
    return out

def main():
    data=json.loads((ROOT/'supports.json').read_text());parts=partitions()
    phase=json.loads((ROOT/'phases.json').read_text())
    assert all(all(d%5 for d in r['diagonal']) for r in phase)
    keep=[r for r in phase if r['binary_feasible']]
    u=variables(BASE,0)
    with (ROOT/'five_input.txt').open('w') as f:
        f.write(str(len(keep))+'\n')
        for r in keep:
            vi,wi=r['vi'],r['wi'];v=variables(data['support_orbits'][vi]['representative'],8);w=variables(parts[wi],16)
            f.write(f'{vi} {wi}\n')
            for a,b,c in zip(u,v,w):f.write(f'{a} {b} {c}\n')
    print('Prepared',len(keep),'quintic label systems.',flush=True)
if __name__=='__main__':main()
