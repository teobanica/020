"""Verify stored identities without rerunning the polynomial GCD or kernel solvers."""
from pathlib import Path
from itertools import combinations,permutations
from concurrent.futures import ProcessPoolExecutor,as_completed
import json,collections,time
from cyclotomic60 import K,Q,POW
from quadratic11 import L
from intervals import real_bound
ROOT=Path(__file__).resolve().parent

def readK(v):return K([Q(x) for x in v])
def readL(v):return L(*[readK(x) for x in v])
def det(a,zero=K()):
 n=len(a);out=zero
 for s in permutations(range(n)):
  v=1
  for i,j in enumerate(s):v=v*a[i][j]
  out=out+((-1)**sum(s[i]>s[j] for i in range(n) for j in range(i+1,n)))*v
 return out

def trim(p):
 while p and not p[-1]:p.pop()
 return p

def add(a,b):return trim([(a[i] if i<len(a) else K())+(b[i] if i<len(b) else K()) for i in range(max(len(a),len(b)))])
def mul(a,b):
 c=[K()] * max(0,len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]=c[i+j]+x*y
 return trim(c)

def decode(line,H):
 w=line.split();blocks=list(map(int,w[1:12]));ex=list(map(int,w[12:]));k=max(blocks)+1
 assert sorted(blocks.count(j) for j in range(k))==sorted({'D':[5,2,2,2],'F':[3,2,2,2,2]}[w[0]])
 assert all(v%7==0 for v in ex)
 for c in range(k):assert sorted(ex[j]%420 for j in range(11) if blocks[j]==c)==list(range(0,420,420//blocks.count(c)))
 coefficients=[[[POW[(ex[j]//7-3*H[r+1][j])%60] for j in range(11) if blocks[j]==c] for c in range(k)] for r in range(4)]
 def matrix(q):
  return [[sum((POW[(ex[j]//7-3*H[r+1][j])%60]*(q if r==3 and j>=5 else 1) for j in range(11) if blocks[j]==c),K()) for c in range(k)] for r in range(4)]
 return k,blocks,ex,matrix

def evaluated_polynomial(k,matrix,pi):
 descriptions=[]
 if k<=4:
  for rs in combinations(range(4),k):descriptions.append(('rank',rs,0))
 for rs in combinations(range(4),k-1):
  for j in range(1,k):descriptions.append(('norm',rs,j))
 kind,rs,j=descriptions[pi]
 if kind=='rank':
  A=matrix(K(0));B=matrix(K(1));a=det([[A[r][c] for c in range(k)] for r in rs]);b=det([[B[r][c] for c in range(k)] for r in rs]);return trim([a,b-a])
 vals=[]
 for q in (K(1),K(-1),POW[15]):
  A=matrix(q);d=[]
  for omit in (0,j):
   z=det([[A[r][c] for c in range(k) if c!=omit] for r in rs]);d.append(z*z.conj())
  vals.append(q*(d[1]-d[0]))
 one,minus,imag=vals;beta=(one-minus)/2;s=(one+minus)/2;d=imag-POW[15]*beta
 return trim([(s+d)/2,beta,(s-d)/2])

def verify_grid():
 total=0;counts=collections.Counter()
 for p in sorted((ROOT/'grid').iterdir()):
  if not (p/'exact_remaining.json').exists():continue
  H=json.loads((p/'seed.json').read_text())['exponents'];lines=(p/'unresolved.txt').read_text().splitlines();records=json.loads((p/'exact_remaining.json').read_text())
  assert len(lines)==len(records)
  for idx,(line,rec) in enumerate(zip(lines,records)):
   k,bl,ex,mat=decode(line,H);A=mat(K(1));N=[[readK(x) for x in row] for row in rec['kernel_basis']];free=rec['free_columns'];assert len(free)==2 and k==5
   assert all(N[f][j]==int(i==j) for i,f in enumerate(free) for j in range(2))
   assert all(sum((A[r][j]*N[j][s] for j in range(k)),K())==0 for r in range(4) for s in range(2))
   assert any(det([[A[r][c] for c in cs] for r in rs]) for rs in combinations(range(4),3) for cs in combinations(range(5),3))
   if rec['exclusion']=='forced_zero_cycle_phase':assert not any(N[rec['column']])
   else:
    assert rec['exclusion']=='inconsistent_linear_norm_equations'
    eq=[[a.conj()*b,a*b.conj(),1-a*a.conj()-b*b.conj()] for a,b in N];w=[readK(v) for v in rec['inconsistency_combination']]
    v=[sum((w[i]*eq[i][j] for i in range(5)),K()) for j in range(3)];assert not v[0] and not v[1] and v[2]
   total+=1;counts[rec['exclusion']]+=1
 return total,dict(counts)

def verify_circle(meta):
 path=ROOT/'circles'/meta['key'];H=json.loads((path/'seed.json').read_text())['exponents'];lines=(path/'unresolved.txt').read_text().splitlines();records=json.loads((path/'exact_certificates.json').read_text());assert len(lines)==len(records)==meta['survivors']
 counts=collections.Counter();pending=set()
 for idx,rec in enumerate(records):
  assert rec['index']==idx;k,bl,ex,matrix=decode(lines[idx],H);g=[readK(x) for x in rec['gcd']];identity=[]
  for term in rec['bezout']:
   f=evaluated_polynomial(k,matrix,term['polynomial_index']);c=[readK(x) for x in term['multiplier']];identity=add(identity,mul(c,f))
  assert identity==g
  case=rec['case'];counts[case]+=1
  if case=='coprime':assert len(g)==1 and g[0]
  elif case=='linear_nonunit':assert len(g)==2 and g[1]==1 and g[0]*g[0].conj()!=1
  elif case in ('no_unit_root','unit_roots_possible'):
   assert len(g)==3 and g[2]==1;c,b=g[:2];assert c*c.conj()==1 and b==c*b.conj();bound=real_bound((b*b.conj()-4).p,rec['digits']);assert bound.serial()==rec['bound']
   if case=='no_unit_root':assert bound.a>0
   else:assert bound.b<0;pending.add(idx)
  elif case=='unit_point':
   e=rec['root60'];assert e in (0,15,30,45);q=POW[e];assert g==[q*q,-2*q,K(1)] or g==[-q,K(1)]
   shifts=meta['shifts'][:];shifts[meta['variable']]=(-e//3)%20
   target=ROOT/'grid'/('_'.join(map(str,shifts)))
   assert lines[idx] in (target/'unresolved.txt').read_text().splitlines()
  else:raise AssertionError(case)
 algebraic=0
 if pending:
  quad=json.loads((path/'quadratic_certificates.json').read_text());assert {(r['index'],r['sign']) for r in quad}=={(i,s) for i in pending for s in (-1,1)}
  for rec in quad:
   i=rec['index'];q=readL(rec['q']);assert q*q.conj()==1;g=[readK(v) for v in records[i]['gcd']];assert q*q+q*g[1]+g[0]==0
   k,bl,ex,mat=decode(lines[i],H);assert k==4
   A=[[sum((L(POW[(ex[j]//7-3*H[r+1][j])%60])*(q if r==3 and j>=5 else 1) for j in range(11) if bl[j]==c),L()) for c in range(k)] for r in range(4)]
   u=[readL(x) for x in rec['kernel_vector']];assert all(x*x.conj()==1 for x in u)
   assert all(sum((A[r][c]*u[c] for c in range(k)),L())==0 for r in range(4))
   assert any(det([[A[r][c] for c in cs] for r in rs],L()) for rs in combinations(range(4),3) for cs in combinations(range(4),3))
   x=[u[bl[j]]*POW[ex[j]//7] for j in range(11)];fac=x[0].inverse();x=[v*fac for v in x];assert x==[readL(v) for v in rec['row']]
   quotient=[x[j]*(L(POW[-3*H[4][j]%60])*(q if j>=5 else 1)) for j in range(11)]
   assert sum(quotient,L())==0
   assert all(not v for z in quotient for a in (z.a,z.b) for j,v in enumerate(a.p) if j not in (0,15))
   assert all(quotient[j]!=-quotient[k] for j in range(11) for k in range(j))
   algebraic+=1
 return {'key':meta['key'],'counts':dict(counts),'quadratic_candidates_checked':algebraic}
if __name__=='__main__':
 from verify_printed_example import verify
 assert verify()
 start=time.time();assert pow(11,30,61)==60
 total,gridcounts=verify_grid();print('grid certificates',total,gridcounts,flush=True)
 summary=json.loads((ROOT/'circles_summary.json').read_text());out=[]
 with ProcessPoolExecutor(max_workers=8) as pool:
  for fut in as_completed([pool.submit(verify_circle,r) for r in summary]):
   out.append(fut.result())
   if len(out)%32==0:print('verified',len(out),'circles; seconds',round(time.time()-start,1),flush=True)
 report={'grid_certificates':total,'grid_exclusions':gridcounts,'circles':len(out),'counts':dict(sum((collections.Counter(r['counts']) for r in out),collections.Counter())),'quadratic_candidates_checked':sum(r['quadratic_candidates_checked'] for r in out),'independent_polynomial_evaluations':True,'bezout_identities':True,'kernel_and_nonregularity_identities':True,'floating_point_mathematics':False,'printed_example_verified':True}
 (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
