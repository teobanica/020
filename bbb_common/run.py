#!/usr/bin/env python3
"""Reproduce the BBB common-support breadth test, or verify saved certificates."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import argparse
import json
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--verify-only', action='store_true')
    p.add_argument('--workers', type=int, default=1, choices=(1, 2))
    args = p.parse_args()
    if sys.flags.optimize:
        p.error('Run without -O or PYTHONOPTIMIZE.')
    if not args.verify_only:
        for mode in ('primary', 'audit'):
            command = ['g++', '-O3', '-std=c++17', '-DGENERIC']
            if mode == 'audit':
                command.append('-DAUDIT')
            subprocess.run(command + [mode + '.cpp', '-o', mode], cwd=ROOT, check=True)

        def scan(mode, profile):
            start = time.monotonic()
            name = f'{mode}_{profile}'
            with (ROOT / (name + '.log')).open('w') as log:
                command = [str(ROOT / mode), 'base_pair.txt', 'head_rows.txt',
                           profile, name + '.txt']
                subprocess.run(command, cwd=ROOT, stdout=log, stderr=log, check=True)
            result = {'mode': mode, 'profile': profile, 'exit': 0,
                      'seconds': round(time.monotonic() - start, 3)}
            print(result, flush=True)
            return result

        results = []
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(scan, mode, profile)
                       for profile in 'BCDEF' for mode in ('primary', 'audit')]
            for future in as_completed(futures):
                results.append(future.result())
        (ROOT / 'execution.json').write_text(json.dumps(results, indent=2) + '\n')
        # A separate first-witness search, explicitly not an exhaustive census.
        subprocess.run(['g++', '-O3', '-std=c++17', '-DPAIR_ONLY', '-DFIRST_WITNESS',
                        'primary.cpp', '-o', 'witness_search'], cwd=ROOT, check=True)
        with (ROOT / 'witness_D.log').open('w') as log:
            subprocess.run([str(ROOT / 'witness_search'), 'base_pair.txt', 'head_rows.txt',
                            'D', 'witness_D.txt'], cwd=ROOT, stdout=log, stderr=log, check=True)
        (ROOT / 'witness_D.txt.writing').replace(ROOT / 'witness_D.txt')
    subprocess.run([sys.executable, 'verify.py'], cwd=ROOT, check=True)


if __name__ == '__main__':
    main()
