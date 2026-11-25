from pathlib import Path
from collections import Counter
from functools import lru_cache
from math import gcd,lcm,factorial
import json,hashlib
from verify_rows import cover,minus,unpack,readrow
R=Path(__file__).resolve().parent
F=R/'fixed_CCC'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def same(a,b):
 assert digest(a)==digest(b),(a,b)
 return digest(a)
def audit_same(name):
 h=json.loads((R/'audit_hashes.json').read_text())[name]
 assert h['primary_sha256']==h['audit_sha256']==digest(R/name)
 return h['primary_sha256']
manifest=json.loads((R/'files_sha256.json').read_text())
for name,h in manifest.items():assert digest(R/name)==h,name
batches=json.loads((F/'batch_audit.json').read_text())
for p,N in [('B',712800),('C',443520),('D',166320),('E',123200),('F',34650)]:
 last=0
 for b in [b for b in batches if b['profile']==p]:
  assert b['first']==last+1 and b['last']>=b['first']
  last=b['last'];assert b['primary_sha256']==b['audit_sha256']
  for mode in ('primary','audit'):
   z=b[mode+'_summary'].split()
   assert int(z[z.index('tests')+1])==last
   assert int(z[z.index('unique')+1])==b['records']
 assert last==N
@lru_cache(maxsize=200000)
def root_cover(v):
 if not v:return True
 a=v[0]
 for p in (2,3,5,7):
  if len(v)<p:continue
  rem=list(v)
  for k in range(p):
   t=(a+840//p*k)%840
   if t not in rem:break
   rem.remove(t)
  else:
   if root_cover(tuple(rem)):return True
 return False
def fixed_check(path,skip,olds):
 dims=Counter();n=0;seen=set()
 for line in path.open():
  z=line.split();g=list(map(int,z[skip:skip+11]));v=list(map(int,z[skip+11:]));assert len(v)==11
  assert tuple(g+v)not in seen;seen.add(tuple(g+v))
  assert set(g)==set(range(max(g)+1));assert int(z[skip-1])==max(g)+1
  for group in set(g):
   jj=[j for j in range(11)if g[j]==group];assert v[jj[0]]==0
   for o in olds:assert root_cover(tuple(sorted((v[j]-o[j])%840 for j in jj))),(path,n,group)
  dims[max(g)]+=1;n+=1
 print(path.name,n,dict(dims),flush=True)
 return {'records':n,'dimensions':dict(sorted(dims.items()))}
seed=json.loads((R/'new_three_parameter_five.json').read_text());E=seed['exponents'];P=seed['powers']
assert len(seed['cycles'])==10
for cert in seed['cycles']:
 i,k=cert['rows'];assert sorted(j for b in cert['cycles']for j in b)==list(range(11))
 for b in cert['cycles']:
  p=len(b);assert p in(2,3,5,7);anchor=b[0]
  for h,j in enumerate(b):
   assert [(P[i][j][t]-P[k][j][t])-(P[i][anchor][t]-P[k][anchor][t])for t in range(3)]==[0]*3
   assert (E[i][j]-E[k][j]-E[i][anchor]+E[k][anchor]-60//p*h)%60==0
assert lcm(60//gcd(60,E[1][5]),60//gcd(60,E[3][7]))==60
old=readrow(R/'family_3286.txt')
for i in range(3):
 for j in range(11):assert old[i][j]==tuple(P[i][j]+[(14*E[i][j])%840])
core=[[14*x for x in E[i]]for i in(0,3,4)]
z=list(map(int,(F/'groups.txt').read_text().split()));assert z[:3]==[1,0,1]
assert [7*x for x in z[3:14]]==core[1] and [7*x for x in z[15:26]]==core[2]
result={'seed_cycles':10,'seed_minimal_level':60}
result['two_core']=fixed_check(F/'pair_families.txt',3,core[:2])
assert result['two_core']['records']==859958
result['three_core']=fixed_check(F/'extensions.txt',2,core)
assert result['three_core']['records']==100788
result['core_audit_sha256']=audit_same('fixed_CCC/extensions.txt')
result['four_row_audit_sha256']=audit_same('five_compatible.txt')
olds=[[(0,0,0,v)for v in row]for row in core]+[old[1]]
n=0;counts=Counter();seen=set()
for line in (R/'five_compatible.txt').open():
 z=list(map(int,line.split()));assert len(z)==67;r,g,v=unpack(z)
 assert tuple(z)not in seen;seen.add(tuple(z))
 for i in range(3):
  assert r[i][i]>=0 and 0<=r[i][3]<840
  if not r[i][i]:assert r[i]==[0]*4
  else:
   assert all(r[i][j]==0 for j in range(i))
   for j in range(i+1,3):
    if r[j][j]:assert 0<=r[i][j]<r[j][j]
 for o in olds:
  for group in set(g):assert cover([minus(v[j],o[j])for j in range(11)if g[j]==group],r),(n,group)
 rank=sum(bool(r[i][i])for i in range(3));counts[(rank,max(g))]+=1;n+=1
 if n%20000==0:print('four-row symbolic verification',n,flush=True)
assert n==118030
result['four_old_rows']={'records':n,'constraint_rank_and_free_row_phases':{str(k):v for k,v in sorted(counts.items())}}
assert (R/'sixth_families.txt').stat().st_size==0
audit_same('sixth_families.txt')
result['sixth_rows']=0
result['hashes']={str(p.relative_to(R)):digest(p)for p in [F/'pair_families.txt',F/'extensions.txt',R/'five_compatible.txt',R/'sixth_families.txt']}
expected=json.loads((R/'verification.json').read_text())
assert json.loads(json.dumps(result))==expected
print(json.dumps(result,indent=2),flush=True)
import runpy
runpy.run_path(str(R/'verify_separation.py'))
print('All delivered identities, counts and audit digests verified.',flush=True)
