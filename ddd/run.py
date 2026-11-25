"""Reproduce the DDD reduction and the exact independent checks."""
from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
def run(*args):subprocess.run(args,cwd=root,check=True)
run(sys.executable,'supports.py')
run('g++','-O3','-std=c++17','phases.cpp','-o','phases')
run('./phases','phase_input.txt','phases.json')
run(sys.executable,'prepare_five.py')
run('g++','-O3','-std=c++17','five.cpp','-o','five')
run('./five','five_input.txt','five.json')
run(sys.executable,'verify.py')
run(sys.executable,'prepare_classes.py')
run('g++','-O3','-std=c++17','classes.cpp','-o','classes')
run('./classes','classes_input.txt','support_classes.json')
run(sys.executable,'verify_classes.py')
run(sys.executable,'check_simple_family.py')
