from pathlib import Path
from collections import Counter,defaultdict
from fractions import Fraction as Q
from math import gcd,lcm
import json,sys
R=Path(__file__).resolve().parent;pilot='--pilot'in sys.argv;dest=R/('pilot_analysis'if pilot else'analysis');dest.mkdir(exist_ok=True)
records=set();complete=[]
for p in 'BCDEF':
 f=R/f'bc_{p}.txt'
 if f.exists():complete.append(p)
 elif pilot:f=f.with_suffix('.txt.writing')
 else:raise ValueError('Missing completed profile '+p)
 if not f.exists():continue
 for l in f.open():
  z=tuple(map(int,l.split()))
  if len(z)!=67:
   assert pilot;continue
  records.add(z)

N=60;M=840

def unpack(z):return z[:12],z[12::5],[z[13+5*j:17+5*j]for j in range(11)]
def normalize(g,e):
 anchors={};labels={};G=[];E=[]
 for group,v in zip(g,e):
  if group not in labels:labels[group]=len(labels);anchors[group]=v
  G.append(labels[group]);E.append((v-anchors[group])%M)
 return tuple(G+E)

general=[];curves=defaultdict(list);points=defaultdict(list);hist=Counter();diag=Counter()
for z in sorted(records):
 H,g,v=unpack(z);rank=sum(bool(H[i])for i in(0,5,10));hist[rank,max(g)]+=1
 assert H[10]==0 and all(x[2]==0 for x in v)
 if rank==0:
  assert max(g)==0;general.append(z)
 elif rank==1:
  a,b,_,e=H[:4]if H[0]else H[4:8];d=gcd(abs(a),abs(b));normal=(a//d,b//d)
  assert normal in ((1,0),(0,1),(1,-1));assert max(g)==0
  for k in range(d):
   q=Q(e+840*k,840*d)%1;assert N%q.denominator==0,(normal,q)
   curves[normal,int(q*N)].append(z)
 else:
  d,h,e,f,r=H[0],H[1],H[3],H[5],H[7];diag[d,f]+=1
  assert max(g)<=1
  for k in range(f):
   beta=Q(r+M*k,M*f)%1
   for j in range(d):
    alpha=((Q(e,M)-h*beta+j)/d)%1
    assert N%alpha.denominator==0 and N%beta.denominator==0,(alpha,beta)
    points[int(alpha*N),int(beta*N)].append(z)

pointsets={};point_counts=Counter()
with(dest/'point_families.txt').open('w')as out:
 for a in range(N):
  for b in range(N):
   applicable=general+curves.get(((1,0),a),[])+curves.get(((0,1),b),[])+curves.get(((1,-1),(a-b)%N),[])+points.get((a,b),[])
   rows=set()
   for z in applicable:
    H,g,v=unpack(z);e=[(x[3]+x[0]*(M//N)*a+x[1]*(M//N)*b)%M for x in v];rows.add(normalize(g,e))
   pointsets[a,b]=rows;point_counts[len(rows)]+=1
   for row in sorted(rows):out.write(' '.join(map(str,[N*a+b,max(row[:11])+1]+list(row)))+'\n')

curve_meta=[];maxdegree=0
curve_keys=[(normal,q)for normal in ((0,1),(1,-1),(1,0))for q in range(N)]
with(dest/'curves.txt').open('w')as out:
 out.write(str(len(curve_keys))+'\n')
 for index,(normal,q)in enumerate(curve_keys):
  extra=curves.get((normal,q),[])
  rows=set()
  for z in general+extra:
   H,g,v=unpack(z);ks=[];es=[]
   for a,b,c,e in v:
    if normal==(1,0):ks.append(b);es.append((e+a*q*(M//N))%M)
    elif normal==(0,1):ks.append(a);es.append((e+b*q*(M//N))%M)
    else:ks.append(a+b);es.append((e+a*q*(M//N))%M)
   assert all(e%14==0 for e in es)
   rows.add(tuple(ks+es))
  rows=sorted(rows);D=max(max(x-y for x,y in zip(a[:11],b[:11]))-min(x-y for x,y in zip(a[:11],b[:11]))for a in rows for b in rows);maxdegree=max(maxdegree,D)
  out.write(f'{index} {normal[0]} {normal[1]} {q} {len(rows)} {D}\n')
  for row in rows:out.write(' '.join(map(str,row))+'\n')
  curve_meta.append({'id':index,'normal':normal,'root_mod60':q,'rows':len(rows),'max_degree':D})
(dest/'generic_records.txt').write_text(''.join(' '.join(map(str,z))+'\n'for z in general))
report={'pilot':pilot,'complete_profiles':complete,'records':len(records),'family_histogram':{str(k):v for k,v in sorted(hist.items())},'diagonals':{str(k):v for k,v in sorted(diag.items())},'generic_records':len(general),'curves':len(curve_keys),'curves_with_additional_records':len(curves),'curve_data':curve_meta,'rank_two_points':len(points),'all_point_graphs':N*N,'point_rows':{str(k):v for k,v in sorted(point_counts.items())},'max_curve_pair_degree':maxdegree}
(dest/'preparation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in ('curve_data','point_rows')},indent=2))
