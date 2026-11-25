from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import json,collections,time
from cyclotomic60 import K,Q,POW
from line_algebra import make,polynomials,trim,divmodp
from intervals import real_bound
ROOT=Path(__file__).resolve().parent
def add(a,b):return trim([(a[i] if i<len(a) else K())+(b[i] if i<len(b) else K()) for i in range(max(len(a),len(b)))])
def mul(a,b):
 c=[K()] * max(0,len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]=c[i+j]+x*y
 return trim(c)
def neg(a):return [-x for x in a]
def xgcd(a,b):
 r0,r1=a,b;s0,s1=[K(1)],[];t0,t1=[],[K(1)]
 while r1:
  q,r=divmodp(r0,r1);r0,r1=r1,r;s0,s1=s1,add(s0,neg(mul(q,s1)));t0,t1=t1,add(t0,neg(mul(q,t1)))
 if not r0:return [],[],[]
 inv=r0[-1].inverse();return [[x*inv for x in p] for p in (r0,s0,t0)]
def serial(p):return [x.serial() for x in p]
def classify(g):
 if len(g)==1:return {'case':'coprime'}
 if len(g)==2:
  q=-g[0]
  if q*q.conj()!=1:return {'case':'linear_nonunit'}
  return {'case':'unit_point','root60':next((i for i,p in enumerate(POW) if p==q),None),'q':q.serial()}
 if len(g)!=3:return {'case':'unresolved_degree','degree':len(g)-1}
 c,b=g[:2]
 if c*c.conj()!=1 or b!=c*b.conj():return {'case':'not_self_inversive'}
 d=b*b.conj()-4
 if not d:
  q=-b/2;assert q*q.conj()==1
  return {'case':'unit_point','root60':next((i for i,p in enumerate(POW) if p==q),None),'q':q.serial()}
 assert d==d.conj()
 for digits in (3,4,6,8,12,16,24):
  interval=real_bound(d.p,digits)
  if interval.a>0 or interval.b<0:break
 else:return {'case':'unresolved_sign'}
 return {'case':'no_unit_root' if interval.a>0 else 'unit_roots_possible','digits':digits,'bound':interval.serial()}
def run(rec):
 p=ROOT/'circles'/rec['key'];H=json.loads((p/'seed.json').read_text())['exponents'];records=[]
 for idx,line in enumerate((p/'unresolved.txt').read_text().splitlines()):
  k,A0,A1=make(line,H);g=[];comb=[];used=[]
  for pi,f in enumerate(polynomials(k,A0,A1)):
   if not f:continue
   g,s,t=xgcd(g,f)
   comb=[mul(s,c) for c in comb]+[t];used.append((pi,f))
   if len(g)==1:break
  active=[{'polynomial_index':pi,'multiplier':serial(c)} for (pi,f),c in zip(used,comb) if c]
  identity=[]
  for (pi,f),c in zip(used,comb):identity=add(identity,mul(f,c))
  assert identity==g
  r={'index':idx,'profile':line[0],'gcd':serial(g),'bezout':active,**classify(g)};records.append(r)
 (p/'exact_certificates.json').write_text(json.dumps(records,separators=(',',':'))+'\n')
 count=dict(collections.Counter(r['case'] for r in records));roots=sorted(set(r['root60'] for r in records if r['case']=='unit_point' and r['root60'] is not None))
 return {'key':rec['key'],'counts':count,'unit_points_root60':roots,'non_fourth_root_cases':[r['index'] for r in records if r['case']=='unit_point' and r['root60'] not in (0,15,30,45)]}
if __name__=='__main__':
 summary=json.loads((ROOT/'circles_summary.json').read_text());out=[];start=time.time()
 with ProcessPoolExecutor(max_workers=8) as pool:
  for fut in as_completed([pool.submit(run,r) for r in summary]):
   r=fut.result();out.append(r)
   if len(out)%16==0:print(len(out),'done',r['counts'],'seconds',round(time.time()-start,1),flush=True)
 (ROOT/'circles_exact_summary.json').write_text(json.dumps(sorted(out,key=lambda x:x['key']),indent=2)+'\n')
 print('TOTALS',dict(sum((collections.Counter(r['counts']) for r in out),collections.Counter())),flush=True)
