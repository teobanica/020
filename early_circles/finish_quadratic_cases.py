from pathlib import Path
from itertools import product
from functools import lru_cache
from math import isqrt
import json,collections
from cyclotomic60 import K,Q,POW
from quadratic11 import L,rrefL
from line_algebra import make
ROOT=Path(__file__).resolve().parent
I=POW[15]
def rsqrt(x):
 a=isqrt(x.numerator);b=isqrt(x.denominator);assert Q(a*a,b*b)==x;return Q(a,b)
def sqrt_gaussian(v):
 assert all(not x for j,x in enumerate(v.p) if j not in (0,15))
 a=v.p[0] if v.p else Q(0);b=v.p[15] if len(v.p)>15 else Q(0)
 n=rsqrt(a*a+b*b)
 if n+a:
  u=rsqrt((n+a)/2);w=b/(2*u)
 else:u=Q(0);w=rsqrt(-a)
 out=K(u)+w*I;assert out*out==v;return out
def cycles(values,p):
 masks=set()
 for origin in values:
  options=[[j for j,v in enumerate(values) if v==origin*POW[(60//p)*k]] for k in range(p)]
  for cols in product(*options):
   if len(set(cols))==p:masks.add(sum(1<<j for j in cols))
 return sorted(masks)
def cover(masks):
 @lru_cache(None)
 def go(rem):
  if not rem:return ()
  bit=rem&-rem
  for mask in masks:
   if mask&bit and mask&rem==mask:
    tail=go(rem^mask)
    if tail is not None:return (mask,)+tail
  return None
 return go((1<<11)-1)
def run_file(path):
 data=json.loads(path.read_text());H=json.loads((path.parent/'seed.json').read_text())['exponents'];lines=(path.parent/'unresolved.txt').read_text().splitlines();out=[]
 for rec in data:
  if rec['case']!='unit_roots_possible':continue
  g=[K([Q(c) for c in x]) for x in rec['gcd']];c,b=g[:2];h=sqrt_gaussian((b*b-4*c)/11)
  idx=rec['index'];k,A0,A1=make(lines[idx],H);assert k==4
  for sign in (-1,1):
   q=L(-b/2,sign*h/2);assert q*q.conj()==1 and q*q+q*b+c==0
   A=[[L(A0[r][j])+(q*A1[r][j] if r==3 else 0) for j in range(k)] for r in range(4)]
   R,piv=rrefL(A);free=[j for j in range(k) if j not in piv];assert len(free)==1
   f=free[0];u=[L(j==f) for j in range(k)]
   for r,j in enumerate(piv):u[j]=-R[r][f]
   assert all(sum((A[r][j]*u[j] for j in range(k)),L())==0 for r in range(4))
   unit=all(v*v.conj()==1 for v in u)
   result={'index':idx,'sign':sign,'q':q.serial(),'rank':len(piv),'free_column':f,'kernel_vector':[v.serial() for v in u],'unit_kernel':unit}
   if unit:
    w=lines[idx].split();blocks=list(map(int,w[1:12]));e=list(map(int,w[12:]));x=[u[blocks[j]]*POW[e[j]//7] for j in range(11)];factor=x[0].inverse();x=[v*factor for v in x]
    old=[[L(POW[3*H[r][j]])*(q.conj() if r==4 and j>=5 else 1) for j in range(11)] for r in range(5)]
    assert all(sum((x[j]*old[r][j].conj() for j in range(11)),L())==0 for r in range(5))
    checks=[]
    for r in range(5):
     values=[x[j]*old[r][j].conj() for j in range(11)];cy={p:cycles(values,p) for p in (2,3,5)};co=cover(sum(cy.values(),[]));checks.append({'row':r,'cycle_masks':cy,'cover':co})
    result['row']=[v.serial() for v in x];result['quotient_checks']=checks;result['regular']=all(a['cover'] is not None for a in checks)
   out.append(result)
 if out:(path.parent/'quadratic_certificates.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
 return out
if __name__=='__main__':
 # 11 is a quadratic nonresidue at a prime that splits in Q(zeta_60).
 prime=next(p for p in (61,181,241,421) if pow(11,(p-1)//2,p)==p-1);print('irreducibility prime',prime,flush=True)
 records=[]
 for path in sorted((ROOT/'circles').glob('*/exact_certificates.json')):
  out=run_file(path)
  if out:
   records.extend(out);print(path.parent.name,len(out),'unit',sum(r['unit_kernel'] for r in out),'regular',sum(r.get('regular',False) for r in out),flush=True)
 print('TOTAL',len(records),'unit',sum(r['unit_kernel'] for r in records),'regular',sum(r.get('regular',False) for r in records),flush=True)
