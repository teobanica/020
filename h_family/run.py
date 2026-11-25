"""Rebuild both exact implementations and exhaust all five profiles."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import subprocess,time,json,sys
ROOT=Path(__file__).resolve().parent
COUNTS={'B':712800,'C':443520,'D':166320,'E':123200,'F':34650}
def build(item):
    name,audit=item
    subprocess.run(['g++','-O3','-std=c++17']+(['-DAUDIT'] if audit else [])+[str(ROOT/(name+'.cpp')),'-o',str(ROOT/(name+('_audit' if audit else '')))],check=True)
def pipeline(item):
    mode,profile=item;folder=ROOT/mode/profile;folder.mkdir(parents=True,exist_ok=True)
    t=time.monotonic();record={'profile':profile}
    for name,args in [('search',[str(COUNTS[profile]),profile]),('refine',[str(folder)]),('bivariate',[str(folder)])]:
        with (folder/(name+'_progress.log')).open('w') as log:
            result=subprocess.run([str(ROOT/(name+('_audit' if mode=='audit' else '')))]+args,cwd=folder,stdout=subprocess.PIPE,stderr=log,text=True)
        (folder/(name+'_summary.txt')).write_text(result.stdout)
        assert result.returncode==0,(mode,profile,name,result.returncode)
        record[name]=result.stdout.strip()
        print(mode,profile,name,result.stdout.strip(),flush=True)
    record['elapsed']=round(time.monotonic()-t,3)
    return record
if __name__=='__main__':
    with ThreadPoolExecutor(6) as pool:list(pool.map(build,[(name,audit)for name in ('search','refine','bivariate')for audit in (False,True)]))
    for mode in ('primary','audit'):
        with ThreadPoolExecutor(5) as pool:records=list(pool.map(pipeline,[(mode,c)for c in COUNTS]))
        (ROOT/(mode+'_summary.json')).write_text(json.dumps(records,indent=2)+'\n')
    subprocess.run([sys.executable,str(ROOT/'verify.py')],check=True)
