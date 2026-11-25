"""Reproduce the complete exceptional-circle calculation using exact arithmetic."""
from pathlib import Path
from itertools import product
from concurrent.futures import ThreadPoolExecutor,as_completed
import subprocess,sys,json,os
ROOT=Path(__file__).resolve().parent
WORKERS=max(1,min(8,os.cpu_count() or 1))
SEED=json.loads((ROOT/'seed.json').read_text())['exponents']
def execute(args,cwd=ROOT):
 p=subprocess.run(args,cwd=cwd,capture_output=True,text=True)
 if p.returncode:raise RuntimeError((args,p.stdout,p.stderr))
 return p.stdout,p.stderr

def grid_case(shifts):
 key='_'.join(map(str,shifts));p=ROOT/'grid'/key;p.mkdir(parents=True,exist_ok=True)
 H=[row[:] for row in SEED]
 for r in range(1,5):
  for j in range(5,11):H[r][j]=(H[r][j]+shifts[r-1])%20
 (p/'seed.json').write_text(json.dumps({'exponents':H}))
 out,err=execute([str(ROOT/'search_grid'),*map(str,shifts)],p)
 assert len(out.splitlines())==6
 (p/'results.txt').write_text(out);(p/'stderr.txt').write_text(err)
 lines=(p/'unresolved.txt').read_text().splitlines()
 assert len(lines)==sum(int(line.split()[5]) for line in out.splitlines()[1:])
 if lines:
  log,_=execute([sys.executable,str(ROOT/'point_exact.py')],p);(p/'exact.log').write_text(log)
 elif (p/'exact_remaining.json').exists():raise AssertionError('Stale point certificate')
 return {'shifts':list(shifts),'survivors':len(lines),'profiles':sorted(set(s[0] for s in lines))}

def circle_case(case):
 var,shifts=case;key=str(var)+'_'+'_'.join(map(str,shifts));p=ROOT/'circles'/key;p.mkdir(parents=True,exist_ok=True)
 H=[SEED[0][:]]
 for r in range(1,5):
  if r!=var+1:H.append([(e+(shifts[r-1] if j>=5 else 0))%20 for j,e in enumerate(SEED[r])])
 H.append(SEED[var+1][:]);(p/'seed.json').write_text(json.dumps({'exponents':H}))
 out,err=execute([str(ROOT/'search_circles'),str(var),*map(str,shifts)],p)
 assert len(out.splitlines())==6
 (p/'results.txt').write_text(out);(p/'stderr.txt').write_text(err)
 lines=(p/'unresolved.txt').read_text().splitlines()
 assert len(lines)==sum(int(line.split()[5]) for line in out.splitlines()[1:])
 return {'key':key,'variable':var,'shifts':shifts,'survivors':len(lines),'profiles':sorted(set(s[0] for s in lines))}

def batch(fun,cases):
 records=[]
 with ThreadPoolExecutor(max_workers=WORKERS) as pool:
  for fut in as_completed([pool.submit(fun,c) for c in cases]):
   records.append(fut.result())
   if len(records)%32==0:print(fun.__name__,len(records),'of',len(cases),flush=True)
 return records

def main():
 for name in ('search_grid','search_circles'):
  execute(['g++','-O3','-std=c++17',name+'.cpp','-o',name])
 grid=batch(grid_case,list(product((0,5,10,15),repeat=4)))
 (ROOT/'grid_summary.json').write_text(json.dumps(sorted(grid,key=lambda r:r['shifts']),indent=2)+'\n')
 cases=[]
 for var in range(4):
  for vals in product((0,5,10,15),repeat=3):
   shifts=list(vals);shifts.insert(var,0);cases.append((var,shifts))
 circles=batch(circle_case,cases)
 (ROOT/'circles_summary.json').write_text(json.dumps(sorted(circles,key=lambda r:r['key']),indent=2)+'\n')
 for name in ('analyze_circles','finish_quadratic_cases','verify_certificates'):
  print('Running',name,flush=True)
  out,err=execute([sys.executable,name+'.py']);(ROOT/(name+'.log')).write_text(out+err)
  print(out.splitlines()[-1],flush=True)
 print('SUCCESS: all 256 coordinate circles excluded for regular sixth rows.',flush=True)
if __name__=='__main__':main()
