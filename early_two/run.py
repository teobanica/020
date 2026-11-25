"""Exhaustive exact reproduction. Python standard library and C++17 only."""
from pathlib import Path
from itertools import combinations,product
from concurrent.futures import ThreadPoolExecutor
import json,subprocess,sys,hashlib
ROOT=Path(__file__).resolve().parent
SEED=json.loads((ROOT/'seed.json').read_text())['exponents']
def cases():
 for var in combinations(range(4),2):
  fixed=[r for r in range(4)if r not in var]
  for vals in product((0,5,10,15),repeat=2):
   shifts=[0]*4
   for r,s in zip(fixed,vals):shifts[r]=s
   yield '_'.join(map(str,(*var,*shifts)))
def compile_all():
 for name in ('search','filter_rank2','isolate_points','filter_cycles'):
  for audit in (False,True):subprocess.run(['g++','-O3','-std=c++17',*(['-DAUDIT']if audit else[]),str(ROOT/(name+'.cpp')),'-o',str(ROOT/(name+('_audit'if audit else'')))],check=True)
def scan(key,audit=False):
 vals=list(map(int,key.split('_')));var=vals[:2];shifts=vals[2:];p=ROOT/('audit'if audit else'surfaces')/key;p.mkdir(parents=True,exist_ok=True)
 H=[SEED[0][:]]
 for r in range(1,5):
  if r-1 not in var:H.append([(e+(shifts[r-1]if j>=5 else 0))%20 for j,e in enumerate(SEED[r])])
 H.extend(SEED[r+1][:]for r in var);(p/'seed.json').write_text(json.dumps({'exponents':H}))
 r=subprocess.run([str(ROOT/('search_audit'if audit else'search')),*key.split('_')],cwd=p,capture_output=True,text=True);assert r.returncode==0,(key,r.stderr)
 (p/'results.txt').write_text(r.stdout);(p/'stderr.txt').write_text(r.stderr)
 rows=r.stdout.splitlines()[1:];assert len(rows)==5
 lines=(p/'unresolved.txt').read_text().splitlines();assert len(lines)==sum(int(l.split()[5])for l in rows)
 if audit:
  orig=ROOT/'surfaces'/key
  assert rows==(orig/'results.txt').read_text().splitlines()[1:],key
  assert (p/'unresolved.txt').read_bytes()==(orig/'unresolved.txt').read_bytes(),key
  for source,dest in [('rank2_unresolved.txt','rank2_unresolved.txt')]:
   if (orig/source).exists():(p/dest).write_bytes((orig/source).read_bytes())
  if (orig/'isolated_points.txt').exists():
   r2=subprocess.run([str(ROOT/'filter_rank2_audit'),*key.split('_')],cwd=p,capture_output=True,text=True);assert not r2.returncode,r2.stderr
   assert r2.stdout==(orig/'rank2_filter.txt').read_text()
   assert (p/'dependencies.txt').read_bytes()==(orig/'dependencies.txt').read_bytes()
   r3=subprocess.run([str(ROOT/'isolate_points_audit'),*key.split('_')],cwd=p,capture_output=True,text=True);assert not r3.returncode,r3.stderr
   assert r3.stdout==(orig/'isolated_points.txt').read_text()
 return {'key':key,'variables':var,'shifts':shifts,'survivors':len(lines),'profiles':sorted(set(l[0]for l in lines))}
def parallel(audit):
 with ThreadPoolExecutor(max_workers=8)as pool:out=list(pool.map(lambda key:scan(key,audit),cases()))
 if not audit:(ROOT/'surfaces_summary.json').write_text(json.dumps(sorted(out,key=lambda x:x['key']),indent=2)+'\n')
 else:(ROOT/'audit_summary.json').write_text(json.dumps({'faces':len(out),'presentations':96*1480490,'agreement':True,'arithmetic':'Gaussian minors, reversed-column Macaulay reduction, Sylvester resultants'},indent=2)+'\n')
 print('audit'if audit else'search','passed',len(out),'faces',flush=True)
if __name__=='__main__':
 if '--audit-only'in sys.argv:parallel(True)
 else:
  compile_all();parallel(False)
  subprocess.run([sys.executable,str(ROOT/'analyze.py')],check=True)
  parallel(True)
  subprocess.run([sys.executable,str(ROOT/'verify_certificates.py')],check=True)
