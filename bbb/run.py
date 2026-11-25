"""Reproduce the BBB support census and its exact checks."""
from pathlib import Path
import subprocess, sys
root=Path(__file__).resolve().parent
def run(*args):subprocess.run(args,cwd=root,check=True)
run(sys.executable,'collapse.py')
run(sys.executable,'prepare_check.py')
run('g++','-O3','-std=c++17','check_collapse.cpp','-o','check_collapse')
run('./check_collapse','check_input.txt','grid_check.json')
run(sys.executable,'prepare_lift.py')
run('g++','-O3','-std=c++17','lift.cpp','-o','lift')
run('./lift','lift_input.txt','lifts.json')
run(sys.executable,'verify.py')
