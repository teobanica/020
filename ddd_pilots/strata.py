"""Stratify a two-dimensional parameter torus into plane, lines and points."""
from collections import defaultdict,Counter
from itertools import combinations
from math import gcd
from pathlib import Path
import json,time
ROOT=Path(__file__).resolve().parent
L=3360
def roots(a,e):
    if a<0:a,e=-a,-e
    assert L%a==0 and e%a==0,(a,e)
    return tuple((e//a+i*(L//a))%L for i in range(a))
def lines(eq):
    a,b,k,e=eq
    assert k==0
    g=gcd(a,b)
    assert a>=0 and (a>0 or b>0)
    return tuple((a//g,b//g,w) for w in roots(g,4*e))
def points(p,q):
    a,b,e=p;c,d,f=q
    if a==c and b==d:return []
    assert a in (0,1) and c in (0,1)
    if not a:p,q=q,p;a,b,e=p;c,d,f=q
    if not c:return [((e-b*f)%L,f)]
    return [((e-b*y)%L,y) for y in roots(d-b,f-e)]
def holds(line,pt):return (sum(x*y for x,y in zip(line[:2],pt))-line[2])%L==0
def build(record_file, outdir):
    targets=defaultdict(lambda: [[],defaultdict(list),defaultdict(list)])
    total=0
    with open(record_file) as src:
        for i,l in enumerate(src):
            row=tuple(map(int,l.split()));assert len(row)==68 and max(row[13::5])<=1
            h=row[0];eq=[row[j:j+4] for j in (1,5,9) if any(row[j:j+4])]
            plane,curves,dots=targets[h];total+=1
            if not eq:plane.append(i)
            elif len(eq)==1:
                for k in lines(eq[0]):curves[k].append(i)
            else:
                assert len(eq)==2 and eq[0][0] in (1,2)
                for a,b,e in lines(eq[0]):
                    for c,d,f in lines(eq[1]):
                        for pt in points((a,b,e),(c,d,f)):dots[pt].append(i)
    report={}
    with open(outdir/'strata.txt','w') as out:
        for h,(plane,curves,dots) in sorted(targets.items()):
            def emit(dim,base,dx,dy,indices):
                print(h,dim,*base,*dx,*dy,len(indices),*indices,file=out)
            emit(2,(0,0),(1,0),(0,1),plane)
            for (a,b,e),indices in sorted(curves.items()):
                if a==1:emit(1,(e,0),(-b,1),(0,0),plane+indices)
                else:
                    assert (a,b)==(0,1)
                    emit(1,(0,e),(1,0),(0,0),plane+indices)
            full=set(dots)
            for p,q in combinations(curves,2):
                if p[:2]!=q[:2]:full.update(points(p,q))
            for pt in sorted(full):
                indices=plane+dots.get(pt,[])
                for k,ii in curves.items():
                    if holds(k,pt):indices+=ii
                emit(0,pt,(0,0),(0,0),indices)
            report[h]={'plane_records':len(plane),'lines':len(curves),'direct_points':len(dots),'all_points':len(full)}
    return report
if __name__=='__main__':
    d=ROOT/'3250';t=time.monotonic();info=build(d/'primary.txt',d)
    info={'families':info,'seconds':time.monotonic()-t}
    print(info);(d/'strata.json').write_text(json.dumps(info,indent=2)+'\n')
