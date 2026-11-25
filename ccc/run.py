"""Reproduce the CCC support/phase reduction and its exact certificates."""
from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
def run(*args):subprocess.run(args,cwd=root,check=True)
run(sys.executable,'supports.py')
run('g++','-O3','-std=c++17','phase_check.cpp','-o','phase_check')
run('./phase_check','phase_input.txt','phase_check.json')
run(sys.executable,'prepare_five.py')
run('g++','-O3','-std=c++17','five.cpp','-o','five')
run('./five','five_input.txt','five.json')
run(sys.executable,'prepare_three.py')
run('g++','-O3','-std=c++17','three.cpp','-o','three')
run('./three','three_input.txt','three.json')
run(sys.executable,'verify.py')
run(sys.executable,'support_classes.py')
