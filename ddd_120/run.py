"""Reproduce the complete level-120 cube-root DDD slice certificate."""
from pathlib import Path
import argparse,json,os,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parent
SCRIPTS=('prepare.py','verify.py','verify_rows.py')
SOURCES=('full.cpp','audit_lowest.cpp','packed_primary.cpp','packed_audit.cpp','graphs.cpp','audit_expansion.cpp')
INPUTS=('ratio_zeta3.txt','reference_R.txt')

def call(args,cwd,log=None):
 print('+',' '.join(map(str,args)),flush=True)
 if log:
  with(cwd/log).open('w')as out:subprocess.run(list(map(str,args)),cwd=cwd,stdout=out,stderr=subprocess.STDOUT,check=True)
 else:subprocess.run(list(map(str,args)),cwd=cwd,check=True)

def compile_cpp(root,source,target,audit=False):
 args=[os.environ.get('CXX','g++'),'-O3','-std=c++17']
 if audit:args.append('-DAUDIT')
 call(args+[root/source,'-o',root/target],root)

def verify(root):
 (root/'build').mkdir(exist_ok=True)
 compile_cpp(root,'audit_expansion.cpp','build/audit_expansion')
 call([root/'build/audit_expansion',root/'analysis/all_families.txt',root/'analysis/point_families.txt',root/'analysis/rows.txt'],root)
 call([sys.executable,'verify.py'],root)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--verify-only',action='store_true')
 p.add_argument('--reuse-scans',action='store_true',help='retain completed scans and regenerate every subsequent calculation')
 args=p.parse_args()
 if args.verify_only:
  if args.reuse_scans:p.error('choose one mode')
  verify(ROOT);return
 dest=ROOT/'regenerated';dest.mkdir(exist_ok=True)
 for name in INPUTS+SCRIPTS+SOURCES:shutil.copy2(ROOT/name,dest/name)
 if args.reuse_scans:
  for profile in'CDEF':
   for prefix in('primary','audit'):
    for suffix in('.txt','.log'):
     name=f'{prefix}_{profile}{suffix}';shutil.copy2(ROOT/name,dest/name)
 else:
  for source,target,audit in(('full.cpp','full',False),('audit_lowest.cpp','audit_lowest',True),('packed_primary.cpp','packed_primary',False),('packed_audit.cpp','packed_audit',True)):
   compile_cpp(dest,source,target,audit)
  # Serial execution limits memory used by the intermediate row lists.
  for prefix in('primary','audit'):
   for profile in'CDEF':
    exe=('packed_primary'if profile=='F'else'full')if prefix=='primary'else('packed_audit'if profile=='F'else'audit_lowest')
    call([dest/exe,dest/'ratio_zeta3.txt',profile,dest/f'{prefix}_{profile}.txt'],dest,f'{prefix}_{profile}.log')
 call([sys.executable,'prepare.py'],dest,'prepare.log')
 for audit in(False,True):
  suffix='_audit'if audit else'';compile_cpp(dest,'graphs.cpp','graphs'+suffix,audit)
  call([dest/('graphs'+suffix),dest/'analysis/point_families.txt',dest/f'analysis/graphs{suffix}.txt',dest/f'analysis/rows{suffix}.txt'],dest,'graphs'+suffix+'.log')
 verify(dest)
 for name in('analysis/all_families.txt','analysis/point_families.txt','analysis/rows.txt','analysis/rows_audit.txt','analysis/graphs.txt','analysis/graphs_audit.txt'):
  assert(ROOT/name).read_bytes()==(dest/name).read_bytes(),name
 for name in('analysis/preparation.json','analysis/symbolic_verification.json','verification.json'):
  assert json.loads((ROOT/name).read_text())==json.loads((dest/name).read_text()),name
 print('All regenerated data agree with the delivered certificate.',flush=True)
if __name__=='__main__':main()
