"""Independent full-group orbit check of the optimized support canonicalization."""
from itertools import permutations,product
from operator import itemgetter
from pathlib import Path
from collections import Counter
import json,time
ROOT=Path(__file__).resolve().parent

def main():
    start=time.time();perms=list(permutations((1,2,3)));gets=[]
    for axis in permutations(range(3)):
        for pa,pb,pc in product(perms,repeat=3):
            maps=[(0,)+p for p in (pa,pb,pc)];index=[]
            for i,j,k in product(range(4),repeat=3):
                old=[None]*3
                for a,x in enumerate((i,j,k)):old[axis[a]]=maps[a][x]
                index.append(16*old[0]+4*old[1]+old[2])
            assert sorted(index)==list(range(64));gets.append(itemgetter(*index))
    assert len(gets)==1296
    raw={}
    with (ROOT/'classes_input.txt').open() as f:
        n=int(next(f))
        for line in f:
            vi,wi,dim,tor,ok,*cols=map(int,line.split());assert len(cols)==11
            a=bytearray(64)
            for c in cols:a[c]+=1
            key=bytes(a)
            if key in raw:
                assert raw[key][:3]==[dim,tor,ok];raw[key][3]+=1
            else:raw[key]=[dim,tor,ok,1]
    assert sum(r[3] for r in raw.values())==n
    data=json.loads((ROOT/'support_classes.json').read_text());classes=data['classes'];dims=Counter();tors=Counter();covered=0;feasible=0
    for r in classes:
        key=bytes(r['tensor']);orb={bytes(get(key)) for get in gets}
        assert min(orb)==key
        count=0
        for t in orb:
            if t in raw:
                a=raw.pop(t)
                assert a[:3]==[r['dimension'],r['torsion_order'],int(r['feasible'])]
                count+=a[3]
        assert count==r['presentations'] and count>0
        covered+=count
        if r['feasible']:feasible+=1;dims[r['dimension']]+=1;tors[r['torsion_order']]+=1
    assert not raw and covered==n and feasible==data['feasible_support_classes']
    report={'verified':True,'group_size':1296,'presentations':n,'total_support_classes':len(classes),
        'feasible_support_classes':feasible,'feasible_by_dimension':dict(sorted(dims.items())),
        'feasible_by_torsion_order':dict(sorted(tors.items())),'seconds':round(time.time()-start,3)}
    (ROOT/'classes_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
