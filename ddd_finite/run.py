"""Reproduce the complete finite-presentation DDD extension calculation."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,hashlib,json,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parent
LABEL_FILES=('labels_counts.txt','families.txt','zero_seeds.txt','exceptional_seeds.txt')
PAIR_FILES=('pair_families.txt','pair_summary.txt','extensions.txt')
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def command(args,log,err=None):
 temp=Path(str(log)+'.writing')
 with temp.open('w')as out:
  if err is None:subprocess.run(args,cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,check=True)
  else:
   with err.open('w')as errors:subprocess.run(args,cwd=ROOT,stdout=out,stderr=errors,check=True)
 temp.replace(log)
def check_pair(d):
 lines=(d/'pair_summary.txt').read_text().splitlines();assert len(lines)==5
 total=0
 for l,expected in zip(lines,(712800,443520,166320,123200,34650)):
  z=l.split();assert int(z[z.index('tests')+1])==expected;total+=int(z[z.index('unique')+1])
 count=0
 for line in(d/'pair_families.txt').open():assert len(line.split())==25;count+=1
 assert count==total
 return count
def check_extensions(d,summary):
 z=(d/summary).read_text().splitlines()[-1].split();assert z[:2]==['summary','group']
 count=0;seeds=set()
 for line in(d/'extensions.txt').open():
  a=line.split();assert len(a)==24;seeds.add(int(a[0]));count+=1
 assert count==int(z[z.index('extension_families')+1]) and len(seeds)==int(z[z.index('extendible')+1])
 return count
def work_group(g):
 primary=ROOT/'zero_results'/str(g);audit=ROOT/'zero_audit'/str(g)
 for mode,d in(('',primary),('_audit',audit)):
  d.mkdir(parents=True,exist_ok=True)
  command([str(ROOT/('pair'+mode)),str(g),str(ROOT/'zero_groups.txt'),str(d/'pair_families.txt')],d/'pair_summary.txt',d/'pair_progress.log');check_pair(d)
  name='extend_summary.txt'if mode else'extension_summary.txt'
  command([str(ROOT/('extend'+mode)),str(g),str(ROOT/'zero_groups.txt'),str(d/'pair_families.txt'),str(d/'extensions.txt')],d/name,d/'extension_progress.log');check_extensions(d,name)
 for name in PAIR_FILES:assert digest(primary/name)==digest(audit/name),(g,name)
 for suffix in('','_audit'):
  command([str(ROOT/('possible_graphs'+suffix)),str(primary/'extensions.txt'),str(primary/('possible_graphs'+suffix+'.txt'))],primary/('possible_graphs'+suffix+'_summary.txt'))
 assert (primary/'possible_graphs.txt').read_bytes()==(primary/'possible_graphs_audit.txt').read_bytes()
 print('finished group',g,flush=True)
 return str(g),{name:digest(primary/name)for name in PAIR_FILES}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);ap.add_argument('--check-existing',action='store_true');args=ap.parse_args()
 if args.check_existing:subprocess.run([sys.executable,str(ROOT/'verify.py')],check=True);return
 assert args.workers>0
 logs=ROOT/'logs';logs.mkdir(exist_ok=True)
 command([sys.executable,'prepare.py'],logs/'prepare.log')
 for name in('labels','pair','extend','possible_graphs','group_all'):
  for suffix in(('',)if name=='group_all'else('','_audit')):
   cmd=['g++','-std=c++17','-O3']+(['-DAUDIT']if suffix else[])+[name+'.cpp','-o',name+suffix]
   command(cmd,logs/(name+suffix+'_build.log'))
 stages=[]
 for suffix in('','_audit'):
  d=ROOT/('regenerated_labels'+suffix);d.mkdir(exist_ok=True);shutil.copy2(ROOT/'labels_input.txt',d/'labels_input.txt')
  command([str(ROOT/('labels'+suffix)),str(d)],d/'labels_summary.txt',d/'labels_progress.log')
  for name,n in zip(LABEL_FILES,(452,201248,54836,11264)):assert sum(1 for _ in(d/name).open())==n
  stages.append(d)
 labelhashes={}
 for name in LABEL_FILES:
  assert digest(stages[0]/name)==digest(stages[1]/name);labelhashes[name]=digest(stages[0]/name);shutil.copy2(stages[0]/name,ROOT/name)
 (ROOT/'audit_label_hashes.json').write_text(json.dumps(labelhashes,indent=2)+'\n')
 shutil.copy2(stages[0]/'labels_summary.txt',ROOT/'labels_summary.txt')
 command([str(ROOT/'group_all'),str(ROOT)],logs/'grouping.log')
 with ThreadPoolExecutor(args.workers)as pool:hashes=dict(pool.map(work_group,range(36)))
 (ROOT/'audit_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
 command([sys.executable,'verify.py'],ROOT/'verify.log')
 print((ROOT/'verification.json').read_text())
if __name__=='__main__':main()
