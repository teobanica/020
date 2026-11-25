"""Verify identities in the supplied certificates without running the kernel solver."""
from pathlib import Path
from itertools import combinations,permutations
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
from fractions import Fraction
import json,subprocess
from cyclotomic60 import K,POW
ROOT=Path(__file__).resolve().parent

def unpack(v):return K([Fraction(x)for x in v])
def determinant(a):
 if len(a)==1:return a[0][0]
 return sum(((-1)**j*a[0][j]*determinant([[row[k]for k in range(len(a))if k!=j]for row in a[1:]])for j in range(len(a))),K())
def matrix(line,H,q=K(1),s=K(1),mask=None):
 z=line.split();b=list(map(int,z[1:12]));e=list(map(int,z[12:]));assert max(b)==4 and all(x%7==0 for x in e)
 out=[[K()for _ in range(5)]for _ in range(4)]
 for r in range(4):
  for j in range(11):
   if mask is not None and r>=2 and bool(mask&(1<<(r-2)))!=(j>=5):continue
   c=POW[(e[j]//7-3*H[r+1][j])%60]
   if mask is None and j>=5:c=c*(q if r==2 else s if r==3 else 1)
   out[r][b[j]]=out[r][b[j]]+c
 return out

def verify(key):
 d=ROOT/'surfaces'/key;H=json.loads((d/'seed.json').read_text())['exponents'];lines=(d/'unresolved.txt').read_text().splitlines()
 # The independent binaries use Gaussian minors and Sylvester determinants.
 a=ROOT/'audit_aux'/key;a.mkdir(parents=True,exist_ok=True);(a/'unresolved.txt').write_bytes((d/'unresolved.txt').read_bytes())
 for exe,name in [('filter_rank2','rank2_filter.txt'),('isolate_points','isolated_points.txt'),('filter_cycles','cycle_filter.txt')]:
  r=subprocess.run([str(ROOT/(exe+'_audit')),*key.split('_')],cwd=a,capture_output=True,text=True);assert r.returncode==0,(key,exe,r.stderr)
  assert r.stdout==(d/name).read_text(),(key,name)
 for name in ['rank2_unresolved.txt','dependencies.txt','curve_unresolved.txt']:assert (a/name).read_bytes()==(d/name).read_bytes()
 depids=set()
 for line in (d/'dependencies.txt').read_text().splitlines():
  vals=list(map(int,line.split()));idx=vals[0];depids.add(idx)
  mats=[matrix(lines[idx],H,mask=m)for m in range(4)]
  for left,right in zip(vals[1::2],vals[2::2]):
   cols=[c for c in range(5)if c not in (left,right)]
   for rows in combinations(range(4),3):
    for mat in mats:assert not determinant([[mat[r][c]for c in cols]for r in rows]),(key,idx)
 pointids=set();roots=set()
 isolated=[tuple(map(int,x.split()))for x in (d/'isolated_points.txt').read_text().splitlines()]
 records=json.loads((d/'point_certificates.json').read_text());assert len(records)==len(isolated)
 for rec,(idx,qe,se)in zip(records,isolated):
  assert (rec['index'],rec['q_exponent420'],rec['s_exponent420'])==(idx,qe,se)
  assert qe%7==se%7==0;roots.update((qe,se));pointids.add(idx)
  A=matrix(lines[idx],H,POW[qe//7],POW[se//7]);N=[[unpack(x)for x in row]for row in rec['kernel_basis']]
  assert len(N)==5 and all(len(row)==2 for row in N)
  assert all(sum((A[r][j]*N[j][t]for j in range(5)),K())==0 for r in range(4)for t in range(2))
  free=rec['free_columns'];assert N[free[0]]==[K(1),K(0)]and N[free[1]]==[K(0),K(1)]
  rows,cols=rec['rank_minor'];assert len(rows)==len(cols)==3 and determinant([[A[r][c]for c in cols]for r in rows])
  assert rec['exclusion']=='forced_zero_cycle_phase'and not any(N[rec['column']])
 curveids=set(map(int,(d/'curve_unresolved.txt').read_text().split()))
 branches=Counter()
 for line in (d/'cycle_filter.txt').read_text().splitlines():
  idx,ok,total,failed,*counts=map(int,line.split());assert idx in curveids and ok==1 and failed==0 and len(counts)==4 and sum(counts)==total
  branches.update(dict(zip((2,3,5,7),counts)))
 assert len((d/'cycle_filter.txt').read_text().splitlines())==len(curveids)
 assert not(depids&pointids or depids&curveids or pointids&curveids)
 assert depids|pointids|curveids==set(range(len(lines))),key
 return {'key':key,'dependencies':len(depids),'points':len(pointids),'curves':len(curveids),'cycle_branches':dict(branches),'parameter_exponents420':sorted(roots)}
if __name__=='__main__':
 for exe in ('filter_rank2','isolate_points','filter_cycles'):
  if not (ROOT/(exe+'_audit')).exists():subprocess.run(['g++','-O3','-std=c++17','-DAUDIT',str(ROOT/(exe+'.cpp')),'-o',str(ROOT/(exe+'_audit'))],check=True)
 keys=[x['key']for x in json.loads((ROOT/'surfaces_summary.json').read_text())]
 assert len(keys)==96 and len(set(keys))==96
 with ProcessPoolExecutor(max_workers=8)as pool:out=list(pool.map(verify,keys))
 report={'faces':len(out),'labelled_presentations':96*1480490,'stage_one_retained':sum(x['dependencies']+x['points']+x['curves']for x in out),'dependent_column_exclusions':sum(x['dependencies']for x in out),'point_exclusions':sum(x['points']for x in out),'cycle_exclusions':sum(x['curves']for x in out),'cycle_branches':dict(sum((Counter(x['cycle_branches'])for x in out),Counter())),'parameter_exponents420':sorted(set(e for x in out for e in x['parameter_exponents420'])),'full_arithmetic_rerun':json.loads((ROOT/'audit_summary.json').read_text()),'per_face':out}
 assert (report['stage_one_retained'],report['dependent_column_exclusions'],report['point_exclusions'],report['cycle_exclusions'])==(4368,3840,512,16)
 (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items()if k!='per_face'},indent=2))
