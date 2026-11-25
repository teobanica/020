from collections import defaultdict
from itertools import product
from pathlib import Path
import json
from normalize import ROOT,U
def prepare():
 D=ROOT/'3297';main=[tuple(map(int,l.split())) for l in (D/'pair.txt').read_text().splitlines()]
 heads=[tuple(map(int,l.split())) for l in (D/'thirds.txt').read_text().splitlines()]
 found={}
 for n in range(40):
  u=[(2*e//21+n*s)%40 for s,p,q,e in main];assert all(e%21==0 for s,p,q,e in main)
  for h in range(6):
   v=[2*e//21 for s,p,q,e in heads[h*11:h*11+11]];g=[p+2*q for s,p,q,e in heads[h*11:h*11+11]]
   best=None
   for sign in (1,-1):
    vv=[(b-(a if sign<0 else 0))%40 for a,b in zip(u,v)]
    for base in set(u):
     fixed=[sign*(a-base)%40 for a in u]
     if sorted(fixed)!=list(U):continue
     for o in product(*(set(vv[j] for j in range(11) if g[j]==k) for k in range(3))):
      cols=sorted((fixed[j],g[j],(vv[j]-o[g[j]])%40) for j in range(11))
      key=tuple(x for _,a,b in cols for x in (a,b))
      if best is None or key<best:best=key
   if best is not None:found[n,h]=best
 keys=sorted(set(found.values()));directory=D/'slices';directory.mkdir(exist_ok=True)
 with (directory/'target.txt').open('w') as f:
  for i,key in enumerate(keys):print(i,*key[1::2],*key[::2],file=f)
 report={'parameter_exponents_mod40':sorted(set(n for n,h in found)), 'normalized_count':len(keys),
 'assignments':[[n,h,keys.index(k)] for (n,h),k in sorted(found.items())]}
 (directory/'mapping.json').write_text(json.dumps(report,indent=2)+'\n');print(report['parameter_exponents_mod40'],len(keys))
if __name__=='__main__':prepare()
