"""Expand the complete modular records at every pair of old parameters."""
from pathlib import Path
from collections import Counter,defaultdict
import json,sys
R=Path(__file__).resolve().parent;M=120
partial='--partial'in sys.argv
D=R/('pilot'if partial else'analysis');D.mkdir(exist_ok=True)
records=set();completed=[]
for p in'CDEF':
 f=R/f'primary_{p}.txt'
 if f.exists():completed.append(p)
 elif partial:f=f.with_suffix('.txt.writing')
 else:raise ValueError('Incomplete profile '+p)
 if not f.exists():continue
 for line in f.read_text().splitlines():
  z=tuple(map(int,line.split()))
  if len(z)==67:records.add(z)
  else:assert partial
roots={(d,e):tuple(a for a in range(M)if(d*a-e)%M==0)for d in range(1,121)if M%d==0 for e in range(M)}
data=[set()for _ in range(M*M)];diagonals=Counter();free=Counter();expanded_incidence=0
for r in sorted(records):
 a,h,c,e=r[:4];u,b,w,f=r[4:8]
 assert a>0 and b>0 and not(c or u or w)and not any(r[8:12])
 assert M%a==M%b==0
 g=r[12::5];maxg=max(g);free[maxg]+=1;diagonals[a,b]+=1;count=0
 for beta in roots[b,f]:
  for alpha in roots[a,(e-h*beta)%M]:
   count+=1;groups={};anchors={};G=[];E=[]
   for j in range(11):
    group,p,q,z,coeff=r[12+5*j:17+5*j];assert z==0
    value=(coeff+p*alpha+q*beta)%M
    if group not in groups:groups[group]=len(groups);anchors[group]=value
    G.append(groups[group]);E.append((value-anchors[group])%M)
   data[M*alpha+beta].add(tuple(G+E))
 assert count==a*b,(r[:12],count)
 expanded_incidence+=count
with(D/'point_families.txt').open('w')as out:
 for point,rows in enumerate(data):
  for row in sorted(rows):out.write(' '.join(map(str,[point,max(row[:11])+1]+list(row)))+'\n')
(D/'all_families.txt').write_text(''.join(' '.join(map(str,z))+'\n'for z in sorted(records)))
report={'partial':partial,'complete_profiles':completed,'modulus':M,'records':len(records),
 'record_diagonals':{str(k):n for k,n in sorted(diagonals.items())},'free_row_parameters':dict(sorted(free.items())),
 'record_point_incidences':expanded_incidence,'point_family_total':sum(map(len,data)),
 'point_family_histogram':dict(sorted(Counter(map(len,data)).items())),
 'points_with_free_rows':{i:sum(max(r[:11])>0 for r in rows)for i,rows in enumerate(data)if any(max(r[:11])>0 for r in rows)}}
(D/'preparation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items()if k!='points_with_free_rows'},indent=2),flush=True)
