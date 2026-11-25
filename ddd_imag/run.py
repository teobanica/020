"""Reproduce the imaginary-surface calculation without replacing its evidence."""
from pathlib import Path
import argparse,json,os,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parent
SCRIPTS=('inspect_records.py','prepare_graphs.py','run_graphs.py','verify.py',
         'verify_rows.py','verify_extra.py','verify_triangles.py')
SOURCES=('full.cpp','audit_lowest.cpp','packed_primary.cpp','packed_audit.cpp')
INPUTS=('b_equals_ic.txt','reference_R.txt')

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
 compile_cpp(root,'graphs/audit_expansion.cpp','build/audit_expansion')
 call([root/'build/audit_expansion',root/'all_families.txt',root/'analysis/point_families.txt',root/'analysis/curves.txt',120],root)
 call([sys.executable,'verify.py'],root)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--verify-only',action='store_true')
 p.add_argument('--reuse-scans',action='store_true',help='retain delivered completed profile scans and regenerate every later stage')
 args=p.parse_args()
 if args.verify_only:
  if args.reuse_scans:p.error('choose one mode')
  verify(ROOT);return
 dest=ROOT/'regenerated';dest.mkdir(exist_ok=True)
 for name in INPUTS+SCRIPTS+SOURCES:shutil.copy2(ROOT/name,dest/name)
 (dest/'graphs').mkdir(exist_ok=True)
 for source in(ROOT/'graphs').glob('*.cpp'):shutil.copy2(source,dest/'graphs'/source.name)
 if args.reuse_scans:
  for profile in'BCDEF':
   for prefix in('primary','audit'):
    for suffix in('.txt','.log'):
     name=f'{prefix}_{profile}{suffix}';shutil.copy2(ROOT/name,dest/name)
 else:
  for source,target,audit in(('full.cpp','full',False),('audit_lowest.cpp','audit_lowest',True),('packed_primary.cpp','packed_primary',False),('packed_audit.cpp','packed_audit',True)):
   compile_cpp(dest,source,target,audit)
  # Serial scans bound peak memory: F holds millions of intermediate keys.
  for prefix in('primary','audit'):
   for profile in'BCDEF':
    exe=('packed_primary'if profile=='F'else'full')if prefix=='primary'else('packed_audit'if profile=='F'else'audit_lowest')
    call([dest/exe,dest/'b_equals_ic.txt',profile,dest/f'{prefix}_{profile}.txt'],dest,f'{prefix}_{profile}.log')
 for script in('inspect_records.py','prepare_graphs.py'):call([sys.executable,script],dest,script[:-3]+'.log')
 for name,source in(('points','possible_graphs.cpp'),('curves','curve_graphs.cpp'),('triangles','triangle_loci.cpp')):
  for audit in(False,True):compile_cpp(dest,'graphs/'+source,'graphs/'+name+('_audit'if audit else''),audit)
 compile_cpp(dest,'graphs/audit_expansion.cpp','graphs/audit_expansion')
 call([sys.executable,'run_graphs.py'],dest,'graphs.log')
 call([sys.executable,'verify.py'],dest,'verify.log')
 compare=['all_families.txt','analysis/generic_records.txt','analysis/curves.txt',
          'analysis/point_families.txt','analysis/points_graphs.txt','analysis/curves_graphs.txt',
          'analysis/triangle_loci.txt','analysis/triangle_contradictions.txt',
          'analysis/points_graphs_audit.txt.witnesses']
 for name in compare:assert (ROOT/name).read_bytes()==(dest/name).read_bytes(),name
 for name in('verification.json','extra_verification.json','loci.json','analysis/preparation.json','incompatible_triangle_example.json'):
  assert json.loads((ROOT/name).read_text())==json.loads((dest/name).read_text()),name
 print('All regenerated data agree with the delivered certificate.',flush=True)
if __name__=='__main__':main()
