from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent
CASE_COUNTS={3203:32,3217:16}
U=(0,0,2,10,10,18,20,20,26,30,34)
def canonical(fields):
    rows=[[(0,0)]*11]+[list(zip(fields[6+11*i:17+11*i],fields[28+11*i:39+11*i])) for i in (0,1)]
    best=None
    for i in range(3):
        for j in range(3):
            if j==i:continue
            k=3-i-j
            pair=[(x-y,d-e) for (x,d),(y,e) in zip(rows[j],rows[i])]
            if len({d for x,d in pair})!=1:continue
            third=[(x-y,d-e) for (x,d),(y,e) in zip(rows[k],rows[i])]
            charge=[d for x,d in third];freq=Counter(charge)
            assert sorted(freq.values())==[2,9] and max(freq)-min(freq)==1
            q=[int(freq[d]==2) for d in charge]
            for base in {x for x,d in pair}:
                u=[(x-base)%40 for x,d in pair]
                if sorted(u)!=list(U):continue
                for e0,_ in [v for h,v in enumerate(third) if not q[h]]:
                    for e1,_ in [v for h,v in enumerate(third) if q[h]]:
                        cols=sorted((u[h],q[h],(third[h][0]-(e1 if q[h] else e0))%40) for h in range(11))
                        key=tuple(x for col in cols for x in col[1:])
                        if best is None or key<best:best=key
    assert best is not None
    return best
def select(directory=ROOT, class_id=3203):
    rows = [tuple(map(int,line.split())) for line in (directory/'labels.txt').read_text().splitlines()]
    assert len(rows)==len(set(rows))==960 and all(r[:2]==(class_id,1) for r in rows)
    keys = sorted({canonical(r) for r in rows})
    assert len(keys)==CASE_COUNTS[class_id]
    with (directory/'target.txt').open('w') as output:
        for i,key in enumerate(keys):
            print(i,*key[1::2],*key[::2],file=output)
    print(len(keys),'normalized affine descriptions')

if __name__=='__main__':select()
