from pathlib import Path
from itertools import permutations
from collections import Counter
import json
ROOT=Path(__file__).resolve().parent
U=(0,0,2,10,10,18,20,20,26,30,34)
def normalize(row):
 d=row[1];es=[(0,)*11,row[6:17],row[17:28]]
 slopes=[[(0,)*d]*11, [row[28+j*d:28+(j+1)*d] for j in range(11)],
         [row[28+(j+11)*d:28+(j+12)*d] for j in range(11)]]
 best=None
 for i,j,k in permutations(range(3)):
  u=tuple((a-b)%40 for a,b in zip(es[j],es[i]))
  du={tuple(a-b for a,b in zip(x,y)) for x,y in zip(slopes[j],slopes[i])}
  if len(du)!=1:continue
  v=tuple((a-b)%40 for a,b in zip(es[k],es[i]))
  dv=[tuple(a-b for a,b in zip(x,y)) for x,y in zip(slopes[k],slopes[i])]
  freq=Counter(dv);assert sorted(freq.values())==[2,2,7]
  base=next(w for w in freq if freq[w]==7)
  for other in permutations(w for w in freq if w!=base):
   w1,w2=([a-b for a,b in zip(w,base)] for w in other)
   assert w1[0]*w2[1]!=w1[1]*w2[0]
   groups=[(base,*other).index(w) for w in dv]
   for rotation in set(u):
    fixed=[(a-rotation)%40 for a in u]
    if sorted(fixed)!=list(U):continue
    from itertools import product
    for offset in product(*(set(v[l] for l in range(11) if groups[l]==g) for g in range(3))):
     cols=sorted((fixed[l],groups[l],(v[l]-offset[groups[l]])%40) for l in range(11))
     key=tuple(t for _,g,e in cols for t in (g,e))
     if best is None or key<best:best=key
 assert best is not None
 return best
def labels(num):
 return [tuple(map(int,l.split())) for l in (ROOT/str(num)/'labels.txt').read_text().splitlines()]
def make():
 rows=labels(3250);ks=sorted(set(map(normalize,rows)))
 with (ROOT/'3250'/'target.txt').open('w') as f:
  for i,k in enumerate(ks):print(i,*k[1::2],*k[::2],file=f)
 print(len(rows),len(ks))
if __name__=='__main__':make()
