from pathlib import Path
from collections import Counter,defaultdict
from fractions import Fraction as Q
from math import gcd,lcm
import json,sys
R=Path(__file__).resolve().parent;partial='--partial'in sys.argv
allrows=set();completed=[]
for p in'BCDEF':
 f=R/f'primary_{p}.txt'
 if f.exists():completed.append(p)
 elif partial:f=f.with_suffix('.txt.writing')
 else:raise ValueError('Profile incomplete: '+p)
 if not f.exists():continue
 for l in f.read_text().splitlines():
  z=tuple(map(int,l.split()))
  if len(z)==67:allrows.add(z)
  else:assert partial
hist=Counter();diags=Counter();curves=defaultdict(set);points=set();generic=[];free=[]
for z in sorted(allrows):
 rank=sum(bool(z[i])for i in(0,5,10));groups=z[12::5];hist[rank,max(groups)]+=1
 assert z[10]==z[2]==z[6]==0 and all(z[i]==0 for i in range(15,67,5))
 if rank==0:generic.append(z)
 elif rank==1:
  a,b,c,e=z[:4]if z[0]else z[4:8];d=gcd(abs(a),abs(b));v=(a//d,b//d)
  for k in range(d):curves[v].add(((Q(e,840)+k)/d)%1)
 else:
  d,h,c,e=z[:4];u,f,w,r=z[4:8];assert u==0;diags[d,f]+=1
  for k in range(f):
   beta=((Q(r,840)+k)/f)%1
   for j in range(d):points.add((((Q(e,840)-h*beta+j)/d)%1,beta))
 if max(groups):free.append(z)
res={'complete_profiles':completed,'partial':len(completed)<5,'records':len(allrows),
 'histogram':{str(k):v for k,v in sorted(hist.items())},'diagonals':{str(k):v for k,v in sorted(diags.items())},
 'curves':[{'normal':v,'count':len(q),'root_orders':dict(Counter(a.denominator for a in q)),'angles':[[a.numerator,a.denominator]for a in sorted(q)]}for v,q in sorted(curves.items())],
 'points':len(points),'point_coordinate_order':lcm(*(a.denominator for pt in points for a in pt))if points else 1,
 'generic_rows':len(generic),'generic_coefficient_order':lcm(*(840//gcd(840,z[i])for z in generic for i in range(16,67,5)))if generic else 1,
 'generic_character_patterns':dict(Counter(str(tuple((z[i],z[i+1])for i in range(13,67,5)))for z in generic))}
prefix='partial_'if res['partial']else''
(R/(prefix+'loci.json')).write_text(json.dumps(res,indent=2)+'\n')
if not res['partial']:
 (R/'all_families.txt').write_text(''.join(' '.join(map(str,z))+'\n'for z in sorted(allrows)))
 (R/'generic_records.txt').write_text(''.join(' '.join(map(str,z))+'\n'for z in generic))
print(json.dumps({k:v for k,v in res.items()if k not in('curves','generic_character_patterns')},indent=2))
print('curves',[(v,len(q),lcm(*(a.denominator for a in q)))for v,q in sorted(curves.items())])
