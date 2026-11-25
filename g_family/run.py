"""Reproduce the complete fixed-core and character-equation calculations.

Default: all stages, into ./regenerated. --stage refine uses the delivered
complete fixed-core list and repeats only the two character refinements.
No output from an incomplete batch is passed to a subsequent calculation.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,subprocess,sys,json,hashlib,shutil
ROOT=Path(__file__).resolve().parent
TOTALS={'B':712800,'C':443520,'D':166320,'E':123200,'F':34650}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def command(args,log):
 with log.open('w')as f:subprocess.run(list(map(str,args)),stdout=f,stderr=subprocess.STDOUT,check=True)
def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--workers',type=int,default=2)
 ap.add_argument('--stage',choices=['all','refine'],default='all')
 ap.add_argument('--workdir',type=Path,default=ROOT/'regenerated')
 args=ap.parse_args();assert args.workers>0
 W=args.workdir.resolve();W.mkdir(parents=True,exist_ok=True)
 B=W/'bin';B.mkdir(exist_ok=True);F=W/'fixed_CCC';F.mkdir(exist_ok=True)
 for mode,flags in [('primary',[]),('audit',['-DAUDIT'])]:
  for name,source in [('pair_chunks','fixed_CCC/pair_chunks.cpp'),('extend','fixed_CCC/extend.cpp'),('refine_five','refine_five.cpp')]:
   subprocess.run(['g++','-std=c++17','-O3',*flags,str(ROOT/source),'-o',str(B/(name+'_'+mode))],check=True)
 subprocess.run(['g++','-std=c++17','-O2',str(ROOT/'check_lattice.cpp'),'-o',str(B/'check_lattice')],check=True)
 command([B/'check_lattice'],W/'lattice_check.txt');assert (W/'lattice_check.txt').read_bytes()==(ROOT/'lattice_check.txt').read_bytes()
 if args.stage=='all':
  C=F/'chunks';C.mkdir(exist_ok=True);jobs=[(p,a,min(a+3999,n))for p,n in TOTALS.items()for a in range(1,n+1,4000)]
  def work(job):
   p,first,last=job;paths=[];counts=[]
   for mode in ('primary','audit'):
    target=C/f'{mode}_{p}_{first}.txt';log=target.with_suffix('.log')
    # An existing file can be reused only if its complete-range summary and
    # record count agree. It is then compared byte-for-byte to the other run.
    valid=False
    if target.exists()and log.exists():
     z=log.read_text().splitlines()[-1].split()
     if 'unique'in z and 'tests'in z:
      valid=(int(z[z.index('tests')+1])==last and int(z[z.index('unique')+1])==sum(1 for _ in target.open()))
    if not valid:command([B/('pair_chunks_'+mode),0,ROOT/'fixed_CCC/groups.txt',target,last,first,p],log)
    z=log.read_text().splitlines()[-1].split();count=sum(1 for _ in target.open())
    assert int(z[z.index('tests')+1])==last and int(z[z.index('unique')+1])==count
    paths.append(target);counts.append(count)
   assert sha(paths[0])==sha(paths[1]);assert counts[0]==counts[1]
   print('complete',p,first,last,'records',counts[0],flush=True)
   return {'profile':p,'first':first,'last':last,'records':counts[0],'sha256':sha(paths[0])}
  with ThreadPoolExecutor(args.workers)as pool:batch=list(pool.map(work,jobs))
  (F/'batch_audit.json').write_text(json.dumps(batch,indent=2)+'\n')
  seen=set();parts=[]
  with(F/'pair_families.txt').open('w')as out:
   for p,n in TOTALS.items():
    before=len(seen)
    for first in range(1,n+1,4000):
     for line in(C/f'primary_{p}_{first}.txt').open():
      z=line.split();assert len(z)==25;k=tuple(map(int,z[3:]));assert z[0]==p
      if k not in seen:seen.add(k);out.write(line)
    parts.append({'profile':p,'new_families':len(seen)-before})
  assert len(seen)==859958
  assert sha(F/'pair_families.txt')==sha(ROOT/'fixed_CCC/pair_families.txt')
  for mode in ('primary','audit'):
   command([B/('extend_'+mode),0,ROOT/'fixed_CCC/groups.txt',F/'pair_families.txt',F/('extensions_'+mode+'.txt')],F/('extensions_'+mode+'_summary.txt'))
   assert sum(1 for _ in(F/('extensions_'+mode+'.txt')).open())==100788
   assert sha(F/('extensions_'+mode+'.txt'))==sha(ROOT/'fixed_CCC/extensions.txt')
  core=F/'extensions_primary.txt'
 else:
  core=ROOT/'fixed_CCC/extensions.txt'
  expected=json.loads((ROOT/'audit_hashes.json').read_text())['fixed_CCC/extensions.txt']['primary_sha256'];assert sha(core)==expected
 for mode in ('primary','audit'):
  command([B/('refine_five_'+mode),ROOT/'family_3286.txt',core,W/('sixth_'+mode+'.txt'),W/('compatible_'+mode+'.txt')],W/('refine_'+mode+'_summary.txt'))
  assert (W/('sixth_'+mode+'.txt')).stat().st_size==0
  assert sum(1 for _ in(W/('compatible_'+mode+'.txt')).open())==118030
  assert sha(W/('compatible_'+mode+'.txt'))==sha(ROOT/'five_compatible.txt')
  summary=(W/('refine_'+mode+'_summary.txt')).read_text().splitlines()[-1]
  assert summary.startswith('candidates 100788 pairs 118030 sixth 0 nodes '),summary
 subprocess.run([sys.executable,str(ROOT/'verify.py')],check=True)
 result={'stage':args.stage,'primary_and_audit_agree':True,'core_rows':100788,'four_old_row_records':118030,'sixth_row_records':0}
 (W/'reproduction.json').write_text(json.dumps(result,indent=2)+'\n')
 print('Reproduced',json.dumps(result),flush=True)
if __name__=='__main__':main()
