#!/usr/bin/env python3
"""Prepare, verify or reproduce the order-11 calculations. No third-party modules."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time

ROOT = Path(__file__).resolve().parent
CATALOG = json.loads((ROOT / 'catalog.json').read_text())
CASES = CATALOG['components']
PY = sys.executable

VERIFY = {
    'a_branch': [[PY, 'verify.py']],
    'a_square': [[PY, 'verify.py']],
    'bbb': [[PY, 'verify.py']],
    'bbb_common': [[PY, 'run.py', '--verify-only']],
    'ccc': [[PY, 'five_verify.py'], [PY, 'verify.py']],
    'ccc_common': [[PY, 'run.py', '--verify-only']],
    'ddd': [[PY, 'five_verify.py'], [PY, 'verify.py'],
            [PY, 'verify_classes.py'], [PY, 'check_simple_family.py']],
    'eee': [['node', 'eee.js']],
    'eee_generic': [[PY, 'run.py', '--verify-only']],
    'fff': [[PY, 'run.py', '--check-existing'], ['node', 'verify_families.js'],
            [PY, 'verify_examples.py']],
    'fff_generic': [[PY, 'run.py', '--verify-only']],
    'mixed': [[PY, 'verify.py']],
    'h_family': [[PY, 'verify.py']],
    'g_family': [[PY, 'verify.py'], [PY, 'verify_separation.py']],
    'ddd_finite': [[PY, 'verify.py']],
    'ddd_one': [[PY, 'run.py', '--verify-only']],
    'ddd_pilots': [[PY, 'run.py', '--verify-only']],
    'ddd_generic': [[PY, 'run.py', '--verify-only']],
    'ddd_real': [[PY, 'run.py', '--verify-only']],
    'ddd_imag': [[PY, 'run.py', '--verify-only'], [PY, 'verify_triangles.py']],
    'ddd_120': [[PY, 'run.py', '--verify-only']],
    'early_sixth': [[PY, 'verify.py']],
    'early_circles': [[PY, 'verify_certificates.py']],
    'early_two': [[PY, 'verify_certificates.py']],
    'early_checks': [[PY, 'n11_cycle_multiset_audit.py'],
                     [PY, 'n11_three_row_first_computations.py']],
}
REPRODUCE = {name: [[PY, 'run.py']] for name in CASES}
REPRODUCE.update({
    'a_branch': [[PY, 'compute.py'], [PY, 'enumerate_labels.py'], [PY, 'verify.py'],
                 ['g++', '-O3', '-std=c++17', 'brute_labels.cpp', '-o', 'brute_labels'],
                 ['./brute_labels']],
    'a_square': [[PY, 'generate_rows.py'],
                 ['g++', '-O3', '-std=c++17', 'complete.cpp', '-o', 'complete'],
                 ['./complete', 'rows.txt', 'completions.json'], [PY, 'verify.py']],
    'eee': [['node', 'eee.js', '--write']],
    'early_checks': VERIFY['early_checks'],
})


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def initialize(work):
    if work == ROOT or ROOT in work.parents or work in ROOT.parents:
        raise ValueError('Choose a work directory separate from the delivered code folder.')
    marker = work / '.n11-work.json'
    if work.exists() and any(work.iterdir()) and not marker.is_file():
        raise ValueError('Work directory is not empty and is not an N11 work directory.')
    work.mkdir(parents=True, exist_ok=True)
    if not marker.exists():
        marker.write_text(json.dumps({'project': 'regular11', 'format': 1}) + '\n')
    elif json.loads(marker.read_text()) != {'project': 'regular11', 'format': 1}:
        raise ValueError('Unrecognized work-directory marker.')


def prepare(name, work):
    if name == 'early_checks':
        dest = work / 'historical'
        dest.mkdir(exist_ok=True)
        for source in (ROOT / name).glob('*.py'):
            shutil.copyfile(source, dest / source.name)
        return dest
    if name == 'ddd_finite':
        prepare('ddd', work)
    record = CASES[name]
    dest = work / 'cases' / record['original']
    dest.mkdir(parents=True, exist_ok=True)
    marker = dest / '.prepared.json'
    expected = hashlib.sha256(json.dumps(record['archives'], sort_keys=True).encode()).hexdigest()
    if marker.exists():
        if json.loads(marker.read_text()).get('archive_sha256') != expected:
            raise ValueError(f'{name}: different dataset already prepared; choose a fresh --work directory.')
    else:
        seen = set()
        for packed in record['archives']:
            archive = ROOT / name / packed['name']
            if digest(archive) != packed['sha256']:
                raise ValueError(f'{name}: dataset checksum mismatch.')
            with tarfile.open(archive, 'r:xz') as tar:
                for member in tar:
                    path = Path(member.name)
                    if (not member.isfile() or path.is_absolute() or '..' in path.parts
                            or member.name not in record['data'] or member.name in seen):
                        raise ValueError(f'{name}: invalid archive member {member.name!r}.')
                    target = dest / path
                    if dest.resolve() not in target.resolve().parents:
                        raise ValueError('Archive destination escapes the work directory.')
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with tar.extractfile(member) as source, target.open('wb') as out:
                        shutil.copyfileobj(source, out)
                    info = record['data'][member.name]
                    if target.stat().st_size != info['bytes'] or digest(target) != info['sha256']:
                        raise ValueError(f'{name}: data checksum mismatch for {member.name}.')
                    seen.add(member.name)
        if seen != set(record['data']):
            raise ValueError(f'{name}: dataset has missing files.')
        marker.write_text(json.dumps({'archive_sha256': expected}) + '\n')
    # Edit the delivered source; each invocation copies its current version.
    for relative in record['sources']:
        target = dest / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name / relative, target)
    if name in ('ddd_one', 'ddd_pilots'):
        # Reuse the fixed pair from the finite DDD calculation. Restore only
        # this member, so verifying one support does not unpack every group.
        member = 'zero_results/21/pair_families.txt'
        expected_hash = CASES['ddd_finite']['data'][member]['sha256']
        target = dest / 'g21_pair_families.txt'
        if not target.exists():
            temporary = target.with_suffix('.tmp')
            with tarfile.open(ROOT / 'ddd_finite' / 'data-pairs-20-21.tar.xz') as archive:
                with archive.extractfile(member) as source, temporary.open('wb') as output:
                    shutil.copyfileobj(source, output)
            if digest(temporary) != expected_hash:
                raise ValueError('DDD pair group 21 failed its recorded checksum.')
            temporary.replace(target)
        elif digest(target) != expected_hash:
            raise ValueError('Prepared DDD pair group 21 differs from the original.')
    if name == 'ddd_pilots':
        shutil.copyfile(ROOT / 'ddd_one' / 'circle_bounds.py', dest / 'circle_bounds.py')
    return dest


def check():
    expected = json.loads((ROOT / 'SHA256SUMS.json').read_text())
    bad = [name for name, value in expected.items()
           if not (ROOT / name).is_file() or digest(ROOT / name) != value]
    if bad:
        raise ValueError('Missing or modified delivered files: ' + ', '.join(bad))
    print(f'All {len(expected)} recorded file checksums agree. This is an integrity check.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['list', 'check', 'prepare', 'verify', 'reproduce'])
    parser.add_argument('components', nargs='*', help='short names from list; prepare/verify default to all')
    parser.add_argument('--work', type=Path, default=ROOT.parent / 'code020-work')
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error('Run without -O or PYTHONOPTIMIZE; the mathematical verifiers use assertions.')
    if args.action == 'list':
        for name, record in CASES.items():
            print(f'{name:16} {record["original"]}')
        print('early_checks     earlier exploratory checks')
        return
    if args.action == 'check':
        check()
        return
    if args.action == 'reproduce' and not args.components:
        parser.error('Select a component for regeneration; see README.txt for scope and costs.')
    names = args.components or list(VERIFY)
    unknown = set(names) - set(VERIFY)
    if unknown:
        parser.error('Unknown components: ' + ', '.join(sorted(unknown)))
    work = args.work.expanduser().resolve()
    initialize(work)
    logs = work / 'logs'
    logs.mkdir(exist_ok=True)
    env = os.environ.copy()
    env.pop('PYTHONOPTIMIZE', None)
    report = []
    try:
        for name in dict.fromkeys(names):
            dest = prepare(name, work)
            if args.action == 'prepare':
                print(f'{name}: {dest}', flush=True)
                continue
            commands = (VERIFY if args.action == 'verify' else REPRODUCE)[name]
            log = logs / f'{args.action}-{name}.log'
            started = time.monotonic()
            item = {'component': name, 'action': args.action, 'status': 'running',
                    'log': str(log), 'commands': commands}
            report.append(item)
            print(f'{args.action} {name}; output: {log}', flush=True)
            with log.open('w') as output:
                for cmd in commands:
                    output.write('+ ' + ' '.join(cmd) + '\n')
                    output.flush()
                    stdin = (dest / 'brute_input.txt').open() if cmd == ['./brute_labels'] else None
                    try:
                        result = subprocess.run(cmd, cwd=dest, stdin=stdin, stdout=output,
                                                stderr=subprocess.STDOUT, env=env)
                    finally:
                        if stdin:
                            stdin.close()
                    if result.returncode:
                        item['status'] = 'failed'
                        raise RuntimeError(f'{name}: command failed ({result.returncode}); see {log}')
            item['status'] = 'passed'
            item['seconds'] = round(time.monotonic() - started, 2)
            print(f'{name}: passed ({item["seconds"]} s)', flush=True)
    finally:
        if report:
            (logs / f'{args.action}-results.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, RuntimeError, OSError, tarfile.TarError) as error:
        print(f'Error: {error}', file=sys.stderr)
        sys.exit(1)
