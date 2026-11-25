#!/usr/bin/env python3
"""Reproduce or verify the DDD support classes with fixed pair 21."""
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse
import hashlib
import json
import subprocess
import sys

from normalize import select, CASE_COUNTS
from verify import check

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compile_solver(temp, source, executable, audit=False):
    flags = ['g++', '-O3', '-std=c++17']
    if audit:
        flags.append('-DAUDIT')
    subprocess.run(flags + [str(ROOT/source), '-o', str(temp/executable)], check=True)


def reproduce(temp, directory, executable, inputs, name, extra=()):
    target = temp/name
    with (temp/(name+'.log')).open('w') as log:
        subprocess.run([str(temp/executable), *map(str,inputs), str(target), *extra],
                       cwd=directory, stdout=log, check=True)
    if digest(target) != digest(directory/name):
        raise ValueError(f'Reproduction differs from the retained {directory/name}')
    return {'file':name, 'sha256':digest(target)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-only', action='store_true')
    parser.add_argument('--class-id', type=int, choices=tuple(CASE_COUNTS))
    options = parser.parse_args()
    if sys.flags.optimize:
        parser.error('Run without Python optimization.')
    classes = [options.class_id] if options.class_id else CASE_COUNTS
    with TemporaryDirectory(dir=ROOT) as work:
        temp = Path(work)
        if not options.verify_only:
            compile_solver(temp, 'scan.cpp', 'scan')
            compile_solver(temp, 'scan.cpp', 'scan_audit', audit=True)
            compile_solver(temp, 'bound.cpp', 'bound')
        compile_solver(temp, 'audit_graph.cpp', 'audit_graph')
        for class_id in classes:
            directory = ROOT if class_id == 3203 else ROOT/str(class_id)
            stages = []
            if not options.verify_only:
                select(directory, class_id)
                inputs = [ROOT/'g21_pair_families.txt',directory/'target.txt']
                for exe,name in [('scan','extensions.txt'),('scan_audit','audit_extensions.txt')]:
                    stages.append(reproduce(temp,directory,exe,inputs,name))
                stages.append(reproduce(temp,directory,'bound',[directory/'extensions.txt'],
                                        'graphs.txt',('-1',)))
            stages.append(reproduce(temp,directory,'audit_graph',[directory/'audit_extensions.txt'],
                                    'audit_graphs.txt'))
            check(directory,class_id)
            if not options.verify_only:
                (directory/'execution.json').write_text(json.dumps(stages,indent=2)+'\n')

if __name__=='__main__':
    main()
