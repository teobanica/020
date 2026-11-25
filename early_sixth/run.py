"""Reproduce every certified calculation. Python standard library and g++ only."""
from pathlib import Path
import subprocess,sys,json
root=Path(__file__).resolve().parent
def run(args,output=None):
    if output:
        with (root/output).open('w') as f:subprocess.run(args,cwd=root,stdout=f,check=True)
    else:subprocess.run(args,cwd=root,check=True)
run(['g++','-O3','-std=c++17','search.cpp','-o','search'])
run([str(root/'search')],'results.txt')
run([sys.executable,'exact_remaining.py'],'exact_remaining.log')
run([sys.executable,'verify.py'])
run([sys.executable,'generic_witnesses.py'],'generic_witnesses.log')
rows=[line.split() for line in (root/'results.txt').read_text().splitlines()[1:]]
assert sum(int(r[2]) for r in rows)==1480490
assert [int(r[5]) for r in rows]==[0,0,0,0,14]
assert len(json.loads((root/'generic_witnesses.json').read_text()))==14
print('SUCCESS: fixed seed has no regular sixth row; every cycle presentation has a nonzero family obstruction.')
