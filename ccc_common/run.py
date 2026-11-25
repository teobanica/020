#!/usr/bin/env python3
"""Exact CCC common-block scan and certificate verification."""
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
    p.add_argument('--scan-only', action='store_true')
    p.add_argument('--workers', type=int, choices=(1, 2), default=1)
    args = p.parse_args()
    if sys.flags.optimize:
        p.error('Run without -O or PYTHONOPTIMIZE.')
    if not args.verify_only:
        for mode in ('primary', 'audit'):
            subprocess.run(['g++', '-O3', '-std=c++17', mode + '.cpp', '-o', mode],
                           cwd=ROOT, check=True)

        def scan(mode, profile):
            start = time.monotonic()
            name = f'{mode}_{profile}'
            with (ROOT / (name + '.log')).open('w') as log:
                subprocess.run([str(ROOT / mode), 'pair.txt', 'heads.txt', profile, name + '.txt'],
                               cwd=ROOT, stdout=log, stderr=log, check=True)
            result = {'mode': mode, 'profile': profile, 'exit': 0,
                      'seconds': round(time.monotonic() - start, 3)}
            print(result, flush=True)
            return result

        results = []
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(scan, mode, profile) for profile in 'BCDEF'
                       for mode in ('primary', 'audit')]
            for future in as_completed(futures):
                results.append(future.result())
        (ROOT / 'execution.json').write_text(json.dumps(results, indent=2) + '\n')
    if not args.scan_only:
        subprocess.run([sys.executable, 'verify.py'], cwd=ROOT, check=True)


if __name__ == '__main__':
    main()
