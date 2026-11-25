from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess,json
R=Path(__file__).resolve().parent
D=R/'analysis'
N=json.loads((D/'preparation.json').read_text())['grid_order']
def job(item):
 kind,audit=item;tag='_audit'if audit else''
 source='point_families.txt'if kind=='points'else'curves.txt'
 args=[R/'graphs'/(kind+tag),D/source,D/(kind+'_graphs'+tag+'.txt')]
 if kind=='curves':args.append(N)
 with(D/(kind+tag+'.log')).open('w')as out:subprocess.run(list(map(str,args)),stdout=out,stderr=out,check=True)
 print((kind,tag),(D/(kind+tag+'.log')).read_text().splitlines()[-1],flush=True)
with ThreadPoolExecutor(max_workers=4)as pool:list(pool.map(job,[(k,a)for k in('points','curves')for a in(False,True)]))
for audit in(False,True):
 tag='_audit'if audit else''
 with(D/('triangles'+tag+'.log')).open('w')as out:
  subprocess.run(list(map(str,[R/'graphs'/('triangles'+tag),D/'point_families.txt',D/('points_graphs'+tag+'.txt'),D/('triangle_loci'+tag+'.txt')])),stdout=out,stderr=out,check=True)
 print((D/('triangles'+tag+'.log')).read_text().strip(),flush=True)
with(D/'expansion_audit.log').open('w')as out:
 subprocess.run(list(map(str,[R/'graphs/audit_expansion',R/'all_families.txt',D/'point_families.txt',D/'curves.txt',N])),stdout=out,stderr=out,check=True)
print((D/'expansion_audit.log').read_text().strip(),flush=True)
