from pathlib import Path
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
import json,subprocess
from cyclotomic60 import K,POW,rref
from extract_polynomials import build,det
ROOT=Path(__file__).resolve().parent

def analyze(key):
 d=ROOT/'surfaces'/key;H=json.loads((d/'seed.json').read_text())['exponents'];lines=(d/'unresolved.txt').read_text().splitlines()
 for exe,out in [('filter_rank2','rank2_filter.txt'),('isolate_points','isolated_points.txt'),('filter_cycles','cycle_filter.txt')]:
  r=subprocess.run([str(ROOT/exe),*key.split('_')],cwd=d,capture_output=True,text=True);assert r.returncode==0,(key,exe,r.stderr)
  (d/out).write_text(r.stdout)
 deps=0
 for line in (d/'dependencies.txt').read_text().splitlines():
  z=list(map(int,line.split()));idx=z[0];k,mats=build(lines[idx],H);assert k==5
  for a,b in zip(z[1::2],z[2::2]):
   cols=[c for c in range(5)if c not in (a,b)]
   for rows in combinations(range(4),3):
    for mat in mats:assert not det([[mat[r][c]for c in cols]for r in rows]),(key,idx,a,b)
   deps+=1
 records=[]
 for line in (d/'isolated_points.txt').read_text().splitlines():
  idx,qe,se=map(int,line.split());assert qe%7==se%7==0
  k,mats=build(lines[idx],H);q=POW[qe//7];s=POW[se//7]
  A=[[K()for c in range(k)]for r in range(4)]
  # Fixed rows occur in every coefficient matrix: select one only.
  for r in range(2):A[r]=mats[0][r]
  for c in range(k):
   A[2][c]=mats[0][2][c]+q*mats[1][2][c]
   A[3][c]=mats[0][3][c]+s*mats[2][3][c]
  R,piv=rref(A);free=[j for j in range(5)if j not in piv];assert len(free)==2,(key,idx,len(piv))
  N=[[K(j==f)for f in free]for j in range(5)]
  for i,j in enumerate(piv):N[j]=[-R[i][f]for f in free]
  assert all(sum((A[r][j]*N[j][t]for j in range(5)),K())==0 for r in range(4)for t in range(2))
  witness=next((rows,cols)for rows in combinations(range(4),3)for cols in combinations(range(5),3)if det([[A[r][c]for c in cols]for r in rows]))
  rec={'index':idx,'q_exponent420':qe,'s_exponent420':se,'rank':3,'rank_minor':witness,'free_columns':free,'kernel_basis':[[x.serial()for x in row]for row in N]}
  zero=next((j for j,row in enumerate(N)if not any(row)),None)
  if zero is not None:rec.update(exclusion='forced_zero_cycle_phase',column=zero)
  else:
   eq=[[a.conj()*b,a*b.conj(),1-a*a.conj()-b*b.conj()]for a,b in N]
   augmented=[row+[K(i==j)for j in range(5)]for i,row in enumerate(eq)]
   rr,pp=rref(augmented,2)
   witness=next((row for row in rr if not any(row[:2])and row[2]),None)
   if witness is not None:
    rec.update(exclusion='inconsistent_linear_norm_equations',inconsistency_combination=[x.serial()for x in witness[3:]],nonzero_constant=witness[2].serial())
   else:raise AssertionError(('unresolved',key,idx,len(pp)))
  records.append(rec)
 (d/'point_certificates.json').write_text(json.dumps(records,indent=2)+'\n')
 return {'key':key,'dependencies_checked':deps,'points':len(records),'curves':len((d/'curve_unresolved.txt').read_text().splitlines()),'exclusions':dict(Counter(r['exclusion']for r in records))}
if __name__=='__main__':
 keys=[x['key']for x in json.loads((ROOT/'surfaces_summary.json').read_text())]
 with ProcessPoolExecutor(max_workers=8)as pool:out=list(pool.map(analyze,keys))
 (ROOT/'exact_summary.json').write_text(json.dumps(out,indent=2)+'\n')
 print('faces',len(out),'dependencies',sum(x['dependencies_checked']for x in out),'points',sum(x['points']for x in out))
 print(sum((Counter(x['exclusions'])for x in out),Counter()))
