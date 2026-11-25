#!/usr/bin/env python3
"""Verify or repeat the two higher-dimensional DDD pilot calculations."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse, hashlib, json, subprocess, sys, time
from verify import census, check_refine, geometry, hard_lines, normal3, verify_generic, verify_all
from normalize import ROOT, make
from slices import prepare
from strata import build

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''):h.update(block)
    return h.hexdigest()

def require_same(actual,saved):
    if digest(actual)!=digest(saved):raise ValueError(f'Reproduction differs for {saved}.')

def compile_one(name,source,work,audit=False):
    flags=['g++','-std=c++17','-O3']+(['-DAUDIT'] if audit else [])
    subprocess.run(flags+[str(ROOT/source),'-o',str(work/name)],check=True)

def run_generic(work,workers):
    d=ROOT/'3297'
    for mode in ('primary','audit'):
        compile_one(mode,'generic.cpp',work,mode=='audit')
    def task(spec):
        mode,profile=spec
        target=work/f'{mode}_{profile}.txt'
        with (work/(target.name+'.log')).open('w') as log:
            subprocess.run([str(work/mode),str(d/'pair.txt'),str(d/'thirds.txt'),
                            profile,str(target)],stdout=log,stderr=log,check=True)
        for suffix in ('','.pairs','.exceptional0','.exceptional1'):
            require_same(Path(str(target)+suffix),d/(target.name+suffix))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(task,[(m,p) for p in 'BCDEF' for m in ('primary','audit')]))

def run_surface(work,directory,full):
    if full:
        for mode in ('primary','audit'):
            target=work/(mode+'.txt')
            with (work/(mode+'.log')).open('w') as log:
                subprocess.run([str(work/('refine_'+mode)),str(ROOT/'g21_pair_families.txt'),
                                str(directory/'target.txt'),str(target)],stdout=log,check=True)
            require_same(target,directory/target.name)
        build(directory/'primary.txt',work)
        require_same(work/'strata.txt',directory/'strata.txt')
    for mode in (('primary','audit') if full else ('audit',)):
        target=work/(mode+'_graphs.txt')
        subprocess.run([str(work/('bound_'+mode)),str(directory/(mode+'.txt')),
                        str(directory/'strata.txt'),str(target)],check=True)
        require_same(target,directory/target.name)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--verify-only',action='store_true')
    p.add_argument('--class-id',type=int,choices=(3250,3297))
    p.add_argument('--workers',type=int,choices=(1,2),default=2)
    args=p.parse_args()
    if sys.flags.optimize:p.error('Run with Python assertions enabled.')
    chosen=[args.class_id] if args.class_id else [3250,3297]
    start=time.monotonic()
    with TemporaryDirectory(dir=ROOT) as temp:
        work=Path(temp)
        for mode in (('audit',) if args.verify_only else ('primary','audit')):
            compile_one('bound_'+mode,'bounds.cpp',work,mode=='audit')
            if not args.verify_only:compile_one('refine_'+mode,'stage_two.cpp',work,mode=='audit')
        if 3250 in chosen:
            if not args.verify_only:make()
            run_surface(work,ROOT/'3250',not args.verify_only)
        if 3297 in chosen:
            if not args.verify_only:
                prepare()
                run_generic(work,args.workers)
            run_surface(work,ROOT/'3297/slices',not args.verify_only)
    if len(chosen)==2:verify_all()
    elif chosen[0]==3250:
        from normalize import normalize
        tgt=[tuple(map(int,l.split())) for l in (ROOT/'3250/target.txt').read_text().splitlines()]
        keys=sorted({normalize(row) for row in census(3250)})
        assert tgt==[(i,*k[1::2],*k[::2]) for i,k in enumerate(keys)]
        check_refine(ROOT/'3250');n,peaks=geometry(ROOT/'3250')
        assert n==10251 and not any(peaks.values())
    else:
        normal3(census(3297));prepare();verify_generic()
        directory=ROOT/'3297/slices';check_refine(directory);n,peaks=geometry(directory)
        assert n==4335;hard_lines(directory,peaks)
    print('Verified',chosen,'in',round(time.monotonic()-start,2),'seconds.',flush=True)

if __name__=='__main__':main()
