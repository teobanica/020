"""Reproduce the full generic-triple enumeration and independent audit."""
from pathlib import Path
import argparse,subprocess,sys,shutil,json,time
ROOT=Path(__file__).resolve().parent

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--verify-only',action='store_true');parser.add_argument('--output',type=Path,default=ROOT/'regenerated');args=parser.parse_args()
 if args.verify_only:
  subprocess.run([sys.executable,str(ROOT/'verify_exact.py')],check=True);return
 dest=args.output.resolve();dest.mkdir(parents=True,exist_ok=True)
 for n in ('family_3286.txt','verify_exact.py'):shutil.copyfile(ROOT/n,dest/n)
 start=time.time();summary=[]
 for tag,src,flags in [('3286','primary.cpp',[]),('audit_3286','audit.cpp',['-DAUDIT'])]:
  binary=dest/(tag+'_solver')
  subprocess.run(['g++','-O3','-std=c++17','-DGENERIC',*flags,str(ROOT/src),'-o',str(binary)],check=True)
  for profile in 'BCDEF':
   out=dest/f'{tag}_{profile}.txt';log=dest/f'{tag}_{profile}.log'
   with log.open('w')as f:subprocess.run([str(binary),str(dest/'family_3286.txt'),profile,str(out)],stdout=f,stderr=f,check=True)
   last=next(l for l in log.read_text().splitlines()if l.startswith('tests '));print(tag,profile,last,flush=True);summary.append({'run':tag,'profile':profile,'summary':last})
 subprocess.run([sys.executable,str(dest/'verify_exact.py')],check=True)
 (dest/'reproduction.json').write_text(json.dumps({'success':True,'elapsed_seconds':time.time()-start,'runs':summary},indent=2)+'\n')
if __name__=='__main__':main()
