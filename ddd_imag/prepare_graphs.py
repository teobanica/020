from pathlib import Path
from collections import Counter,defaultdict
from fractions import Fraction as Q
from math import gcd,lcm
from itertools import combinations
import json,sys
R=Path(__file__).resolve().parent;partial='--partial'in sys.argv
D=R/('pilot_analysis'if partial else'analysis');D.mkdir(exist_ok=True)
records=set();completed=[]
for p in'BCDEF':
 f=R/f'primary_{p}.txt'
 if f.exists():completed.append(p)
 elif partial:f=f.with_suffix('.txt.writing')
 else:raise ValueError('Missing completed profile '+p)
 if not f.exists():continue
 for l in f.read_text().splitlines():
  z=tuple(map(int,l.split()))
  if len(z)==67:records.add(z)
  else:assert partial
generic=[];curves=defaultdict(list);points=defaultdict(list);hist=Counter()
def unpack(z):return z[:12],z[12::5],[z[13+5*j:17+5*j]for j in range(11)]
for z in sorted(records):
 H,g,v=unpack(z);rank=sum(bool(H[i])for i in(0,5,10));hist[rank,max(g)]+=1
 assert H[2]==H[6]==H[10]==0 and all(x[2]==0 for x in v)
 if rank==0:
  assert max(g)==0 and all(x[:3]==(0,0,0)for x in v);generic.append(z)
 elif rank==1:
  a,b,c,e=H[:4]if H[0]else H[4:8];d=gcd(abs(a),abs(b));normal=(a//d,b//d)
  assert normal in((0,1),(1,0),(1,-1),(1,-2),(2,-1))and max(g)==0
  for k in range(d):curves[normal,((Q(e,840)+k)/d)%1].append(z)
 else:
  a,b,c,e=H[:4];u,d,w,f=H[4:8];assert u==0 and max(g)<=1
  for k in range(d):
   beta=((Q(f,840)+k)/d)%1
   for j in range(a):points[((Q(e,840)-b*beta+j)/a)%1,beta].append(z)
# Determine and check the finite grid from actual equations, including
# every intersection of nonparallel rank-one components.
intersections=set()
for ((a,b),q),((c,d),s)in combinations(curves,2):
 det=a*d-b*c
 if not det:continue
 for j in range(abs(det)):
  for k in range(abs(det)):
   intersections.add((((d*(q+j)-b*(s+k))/det)%1,((-c*(q+j)+a*(s+k))/det)%1))
N=lcm(60,*(v.denominator for point in set(points)|intersections for v in point),*(q.denominator for normal,q in curves))
assert 840%N==0,(N,'point coefficient group must be enlarged')
step=840//N
def norm(g,ex):
 labels={};anchors={};G=[];E=[]
 for a,e in zip(g,ex):
  if a not in labels:labels[a]=len(labels);anchors[a]=e
  G.append(labels[a]);E.append((e-anchors[a])%840)
 return tuple(G+E)
pc=Counter();free=Counter();total=0
with(D/'point_families.txt').open('w')as out:
 for a in range(N):
  for b in range(N):
   applicable=list(generic)
   for normal in {v for v,q in curves}:
    q=(Q(normal[0]*a+normal[1]*b,N))%1
    applicable+=curves.get((normal,q),[])
   applicable+=points.get((Q(a,N),Q(b,N)),[])
   rows=set()
   for z in applicable:
    H,g,v=unpack(z);rows.add(norm(g,[(e+step*(a*p+b*q))%840 for p,q,c,e in v]))
   pc[len(rows)]+=1;total+=len(rows)
   for row in sorted(rows):
    out.write(' '.join(map(str,[N*a+b,max(row[:11])+1]+list(row)))+'\n')
    if max(row[:11]):free[N*a+b]+=1
metadata=[];orders=set();maxspread=0;differences=set();curve_total=0
with(D/'curves.txt').open('w')as out:
 out.write(str(len(curves))+'\n')
 for index,((A,B),q)in enumerate(sorted(curves)):
  root=int(q*840);assert q*840==root
  if (A,B)==(0,1):a0=0;b0=root;da=1;db=0
  elif A==1:a0=root;b0=0;da=-B;db=1
  else:assert B==-1;a0=0;b0=-root;da=1;db=A
  assert A*a0+B*b0==root and A*da+B*db==0 and gcd(da,db)==1
  rows=set()
  for z in generic+curves[(A,B),q]:
   H,g,v=unpack(z);ks=tuple(p*da+q*db for p,q,c,e in v);es=tuple((e+p*a0+q*b0)%840 for p,q,c,e in v)
   rows.add(ks+es)
  rows=sorted(rows);Ds=set()
  for u,v in combinations(rows,2):
   chars={a-b for a,b in zip(u[:11],v[:11])}
   Ds|={abs(a-b)for a,b in combinations(chars,2)}
  order=lcm(*(840//gcd(840,e)for row in rows for e in row[11:]));spread=max(Ds,default=0)
  orders.add(order);differences|=Ds;maxspread=max(maxspread,spread);curve_total+=len(rows)
  out.write(f'{index} {A} {B} {root} {a0} {b0} {da} {db} {len(rows)} {spread}\n')
  for row in rows:out.write(' '.join(map(str,row))+'\n')
  metadata.append({'id':index,'normal':[A,B],'root_mod840':root,'rows':len(rows),'spread':spread,'coefficient_order':order,'pair_character_differences':sorted(Ds)})
base=lcm(210,*orders);K=lcm(840,N,*(base*d for d in differences))
(D/'generic_records.txt').write_text(''.join(' '.join(map(str,z))+'\n'for z in generic))
report={'partial':partial,'complete_profiles':completed,'records':len(records),'histogram':{str(k):v for k,v in sorted(hist.items())},'grid_order':N,'point_count':N*N,'point_family_count':total,'point_family_histogram':dict(sorted(pc.items())),'free_family_points':dict(sorted(free.items())),'rank_two_support_points':len(points),'curve_intersections':len(intersections),'generic_rows':len(generic),'curve_count':len(curves),'curve_family_count':curve_total,'curve_coefficient_orders':sorted(orders),'max_pair_spread':maxspread,'mixed_cycle_root_base':base,'curve_phase_modulus':K,'curves':metadata}
(D/'preparation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items()if k not in('curves','free_family_points','point_family_histogram')},indent=2))
