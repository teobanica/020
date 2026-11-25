from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess
r=Path(__file__).resolve().parent;d=r/'analysis'
def run(tag):
 name,audit=tag;exe=name+('_audit'if audit else'');source={'global':'generic_records.txt','curves':'curves.txt','points':'point_families.txt'}[name];dest=name+'_graph'+('s'if name!='global'else'')+('_audit'if audit else'')+'.txt'
 with(d/(exe+'.log')).open('w')as f:subprocess.run([str(r/'graphs'/exe),str(d/source),str(d/dest)],stdout=f,stderr=f,check=True)
 print(exe,(d/(exe+'.log')).read_text().splitlines()[-1],flush=True)
with ThreadPoolExecutor(max_workers=4)as pool:list(pool.map(run,[(n,a)for n in ('global','curves','points')for a in (False,True)]))
with(d/'expansion_audit.log').open('w')as f:subprocess.run([str(r/'graphs/audit_expansion'),str(r/'all_families.txt'),str(d/'point_families.txt'),str(d/'curves.txt')],stdout=f,stderr=f,check=True)
print((d/'expansion_audit.log').read_text(),flush=True)
