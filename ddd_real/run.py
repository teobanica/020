"""Reproduce the exceptional-surface certificate; retain delivered data."""
from pathlib import Path
import argparse, json, os, shutil, subprocess, sys

ROOT=Path(__file__).resolve().parent
SCRIPTS=('organize.py','prepare_graphs.py','publish_generic.py',
         'match_final_five.py','run_graphs.py','verify_rows.py',
         'verify_surface.py','verify_extra.py')
INPUTS=('b_equals_c.txt','reference_R.txt','reference_G.json',
        'exceptional_equations_previous.txt')

def call(args,cwd,log=None):
 print('+',' '.join(map(str,args)),flush=True)
 if log:
  with (cwd/log).open('w')as f:
   subprocess.run(list(map(str,args)),cwd=cwd,stdout=f,stderr=subprocess.STDOUT,check=True)
 else:subprocess.run(list(map(str,args)),cwd=cwd,check=True)

def compile_cpp(root,source,target,audit=False):
 args=[os.environ.get('CXX','g++'),'-O3','-std=c++17']
 if audit:args.append('-DAUDIT')
 call(args+[str(root/source),'-o',str(root/target)],root)

def verify(root):
 # This independent expansion audit tests the saved equations directly.
 (root/'build').mkdir(exist_ok=True)
 compile_cpp(root,'graphs/audit_expansion.cpp','build/audit_expansion')
 call([root/'build/audit_expansion',root/'all_families.txt',
       root/'analysis/point_families.txt',root/'analysis/curves.txt'],root)
 call([sys.executable,'verify_surface.py'],root)
 call([sys.executable,'verify_extra.py'],root)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--verify-only',action='store_true',help='verify delivered records and graphs, without re-enumeration')
 p.add_argument('--reuse-scans',action='store_true',help='copy the delivered completed scans; regenerate all later calculations')
 args=p.parse_args()
 if args.verify_only:
  if args.reuse_scans:p.error('choose one mode')
  verify(ROOT);return
 dest=ROOT/'regenerated';dest.mkdir(exist_ok=True)
 for name in INPUTS+SCRIPTS+('full.cpp','audit_lowest.cpp','packed_audit.cpp'):
  shutil.copy2(ROOT/name,dest/name)
 (dest/'graphs').mkdir(exist_ok=True)
 for source in(ROOT/'graphs').glob('*.cpp'):shutil.copy2(source,dest/'graphs'/source.name)
 if args.reuse_scans:
  for profile in 'BCDEF':
   for prefix in ('bc','audit_bc'):
    for suffix in ('.txt','.log'):
     name=f'{prefix}_{profile}{suffix}';shutil.copy2(ROOT/name,dest/name)
 else:
  compile_cpp(dest,'full.cpp','full')
  compile_cpp(dest,'audit_lowest.cpp','audit_lowest',True)
  compile_cpp(dest,'packed_audit.cpp','packed_audit',True)
  # Serial scans avoid overlapping the large intermediate-record sets.
  for prefix in ('bc','audit_bc'):
   for profile in 'BCDEF':
    exe='full'if prefix=='bc'else('packed_audit'if profile=='F'else'audit_lowest')
    call([dest/exe,dest/'b_equals_c.txt',profile,dest/f'{prefix}_{profile}.txt'],dest,f'{prefix}_{profile}.log')
 records=set()
 for profile in 'BCDEF':
  for line in(dest/f'bc_{profile}.txt').open():
   z=tuple(map(int,line.split()));assert len(z)==67;records.add(z)
 (dest/'all_families.txt').write_text(''.join(' '.join(map(str,z))+'\n'for z in sorted(records)))
 for script in ('organize.py','prepare_graphs.py','publish_generic.py','match_final_five.py'):
  call([sys.executable,script],dest,script[:-3]+'.log')
 for name,source in (('global','global_graph.cpp'),('curves','curve_graphs.cpp'),('points','possible_graphs.cpp')):
  for audit in(False,True):
   compile_cpp(dest,'graphs/'+source,'graphs/'+name+('_audit'if audit else''),audit)
 compile_cpp(dest,'graphs/audit_expansion.cpp','graphs/audit_expansion')
 call([sys.executable,'run_graphs.py'],dest)
 call([sys.executable,'verify_surface.py'],dest)
 call([sys.executable,'verify_extra.py'],dest)
 # Primary/audit agreement is checked above. Check delivered/regenerated
 # data as sets where enumeration order is immaterial, bytes otherwise.
 compare=['all_families.txt','analysis/generic_records.txt',
          'analysis/generic_published.txt','analysis/curves.txt',
          'analysis/point_families.txt','analysis/global_graph.txt',
          'analysis/curves_graphs.txt','analysis/points_graphs.txt']
 for name in compare:assert (ROOT/name).read_bytes()==(dest/name).read_bytes(),name
 for name in('verification.json','extra_verification.json','connected_tori.json','analysis/five_equivalences.json'):
  assert json.loads((ROOT/name).read_text())==json.loads((dest/name).read_text()),name
 print('All regenerated certificate data agree with the delivered data.',flush=True)

if __name__=='__main__':main()
