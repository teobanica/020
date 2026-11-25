"""Verify the label censuses, literal covers, torus strata and completion bounds."""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import hashlib,json,sys
from normalize import U, normalize
from slices import prepare
ROOT=Path(__file__).resolve().parent
PAIR_HASH='0d251854aedf2de521bbbbd3c587d6da247207c9db2e084a12c77ebc60f4a496'
L=3360

def mm(a,b):
    return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]

def det(rows):
    a=[[Fraction(v) for v in r] for r in rows];ans=1
    for j in range(len(a)):
        pivot=next((i for i in range(j,len(a)) if a[i][j]),None)
        if pivot is None:return 0
        if pivot!=j:a[j],a[pivot]=a[pivot],a[j];ans=-ans
        ans*=a[j][j]
        for i in range(j+1,len(a)):
            q=a[i][j]/a[j][j]
            a[i]=[x-q*y for x,y in zip(a[i],a[j])]
    return ans

def census(num):
    d=2 if num==3250 else 3;directory=ROOT/str(num)
    c=json.loads((directory/'class.json').read_text());assert c['class_id']==num and c['dimension']==d
    supports=c['supports'];maps=[{j:(i,k) for i,b in enumerate(part) for k,j in enumerate(b)} for part in supports]
    weights=[];sg=[]
    for j in range(11):
        (p,pp),(q,qp),(r,rp)=[m[j] for m in maps]
        row=[0]*6
        if p:row[p-1]=1
        if q:row[2+q]=-1
        weights.append(row);sg.append((rp if r else 0)-(pp if p else 0)+(qp if q else 0))
    edges=[(j,b[0]) for b in supports[2] for j in b[1:]]
    A=[[x-y for x,y in zip(weights[j],weights[h])] for j,h in edges]
    g=[sg[j]-sg[h] for j,h in edges];S,V=c['U'],c['V'];di=c['diagonal'];r=6-d
    assert A==c['A'] and g==c['g'] and abs(det(S))==abs(det(V))==1
    assert di==[1]*(r-1)+[2] and mm(mm(S,A),V)==[[di[i] if i==j and i<r else 0 for j in range(6)] for i in range(7)]
    ug=[sum(x*y for x,y in zip(row,g)) for row in S];assert not any(a%2 for a in ug[r:])
    labs=[(0,)+q for q in permutations(range(1,5))];rows=set()
    for masks in product(range(24),repeat=3):
        head=[]
        for k,part in zip(masks,supports):
            h=[0]*11
            for j,a in zip(part[0],labs[k]):h[j]=a
            head.append(h)
        rhs=[(head[2][j]-head[0][j]+head[1][j])-(head[2][h]-head[0][h]+head[1][h]) for j,h in edges]
        uf=[sum(x*y for x,y in zip(row,rhs)) for row in S]
        if any(x%5 for x in uf[r:]):continue
        for branch in range(2):
            y=[8*(pow(di[i],-1,5)*uf[i]%5)+(20//di[i])*(ug[i]+(2*branch if di[i]==2 else 0)) for i in range(r)]
            phases=[sum(a*b for a,b in zip(row,y)) for row in V];ee=[];sl=[]
            for a in range(2):
                for j in range(11):
                    b,k=maps[a][j];coord=a*3+b-1
                    ee.append((8*head[a][j]+(20*k+phases[coord] if b else 0))%40)
                    sl.extend(V[coord][r:] if b else [0]*d)
            rows.add((num,d,*masks,branch,*ee,*sl))
    saved=[tuple(map(int,l.split())) for l in (directory/'labels.txt').read_text().splitlines()]
    assert len(saved)==len(rows)==144 and rows==set(saved)
    return rows

def normal3(rows):
    matches=Counter()
    for raw in rows:
        u=raw[6:17];v=raw[17:28];du=[raw[28+3*j:31+3*j] for j in range(11)]
        dv=[raw[61+3*j:64+3*j] for j in range(11)]
        r=[(a-b)%40 for a,b in zip(u,v)];y=[-a%40 for a in v]
        dr=[tuple(a-b for a,b in zip(x,t)) for x,t in zip(du,dv)]
        dy=[tuple(-a for a in x) for x in dv]
        assert all(x==(0,0,0) for x in dr[:5]+dy[:5])
        assert sorted(r[:5])==list(range(0,40,8));cols=sorted(range(5),key=lambda j:r[j])
        slope=next(a for a in (2,3,4) if [y[j] for j in cols]==[8*a*i%40 for i in range(5)])
        assert len(set(dr[5:]))==1 and len(set(dy[5:7]))==len(set(dy[7:]))==1
        assert abs(det([dr[5],dy[5],dy[7]]))==1
        for sign,last in [(1,[5,6,7,10,9,8]),(-1,[5,6,7,10,8,9])]:
            if [(r[j]-r[5])%40 for j in last]!=[0,0,20,20,10,30]:continue
            expected=[0,20,0,10,20,30] if sign==1 else [0,20,0,30,10,20]
            assert [(y[j]-(y[5] if j<7 else y[7]))%40 for j in last]==expected
            matches[(slope,sign)]+=1;break
        else:raise ValueError('Unmapped class-3297 component.')
    assert len(matches)==6 and sum(matches.values())==144
    return matches

def cover(row,mod):
    @lru_cache(None)
    def solve(mask):
        if not mask:return True
        h=(mask&-mask).bit_length()-1
        for p in (2,3,5,7):
            if p>mask.bit_count():continue
            choices=[[j for j in range(11) if mask>>j&1 and row[j]==row[h][:-1]+((row[h][-1]+mod//p*k)%mod,)] for k in range(1,p)]
            for other in product(*choices):
                if solve(mask^((1<<h)+sum(1<<j for j in other))):return True
        return False
    return solve(2047)

def verify_generic():
    directory=ROOT/'3297'
    u=[tuple(map(int,l.split())) for l in (directory/'pair.txt').read_text().splitlines()]
    vs=[tuple(map(int,l.split())) for l in (directory/'thirds.txt').read_text().splitlines()]
    assert u==[(0,0,0,84*j) for j in range(5)]+[(1,0,0,k) for k in (0,0,210,210,105,315)]
    expected=[]
    for a in (2,3,4):
        for sign in (1,-1):
            expected.extend([(0,0,0,84*j*a%420) for j in range(5)]+[(0,1,0,0),(0,1,0,210)]+[(0,0,1,e) for e in ((0,105,210,315) if sign==1 else (0,315,105,210))])
    assert vs==expected
    summary={};exceptions=set();finals=defaultdict(set);one=[(0,)*4]*11
    for profile,count in zip('BCDEF',(712800,443520,166320,123200,34650)):
        records=[]
        for mode in ('primary','audit'):
            name=directory/f'{mode}_{profile}.txt'
            pairs={tuple(map(int,l.split())) for l in Path(str(name)+'.pairs').read_text().splitlines()}
            leaves={tuple(map(int,l.split())) for l in name.read_text().splitlines()}
            log=(directory/f'{mode}_{profile}.log').read_text()
            assert f'tests {count} pairs {len(pairs)} final {len(leaves)} nodes' in log
            records.append((pairs,leaves))
        assert records[0]==records[1]
        for hkey in leaves:
            h=hkey[0];row=tuple(hkey[j:j+5] for j in range(1,56,5));finals[h].add(row)
            for v in (one,u,vs[11*h:11*h+11]):
                q=[(x[0],*(x[i+1]-y[i] for i in range(3)),(x[-1]-y[-1])%420) for x,y in zip(row,v)]
                assert cover(q,420)
        for raw in pairs:
            row=[raw[j:j+5] for j in range(0,55,5)]
            for v in (one,u):
                q=[(x[0],*(x[i+1]-y[i] for i in range(3)),(x[-1]-y[-1])%420) for x,y in zip(row,v)]
                assert cover(q,420)
        for stage in (0,1):
            exceptions.update(tuple(map(int,l.split())) for l in (directory/f'audit_{profile}.txt.exceptional{stage}').read_text().splitlines())
        summary[profile]=(count,len(pairs),len(leaves))
    assert len(exceptions)==1982 and all(any(v[:3]) and 0<=v[-1]<420 and max(map(abs,v[:3]))<=5 for v in exceptions)
    (directory/'exceptions.txt').write_text(''.join(' '.join(map(str,v))+'\n' for v in sorted(exceptions)))
    witness={}
    for h,rows in finals.items():
        cs=sorted(r for r in rows if max(x[0] for x in r))
        assert len(rows)==526 and len(cs)==6
        for c in cs:
            assert [x[0] for x in c]==[0]*5+[1]*6 and all(x[-1]%105==0 for x in c[5:])
            assert any(all(x==(0,0,0,0,84*k*j%420) for j,x in enumerate(c[:5])) for k in set(range(5))-{0,1,2+h//2})
        for r in rows:
            if r in cs:continue
            assert any(all(len({tuple(x[i]-y[i] for i in range(1,4))+((x[-1]-y[-1])%420,) for x,y in zip(r,c) if y[0]==g})==1 for g in (0,1)) for c in cs)
        success=None
        for i,j in combinations(range(6),2):
            x,y=cs[i],cs[j]
            q=[tuple(a[k]-b[k] for k in range(1,4))+((a[-1]-b[-1])%420,) for a,b in zip(x,y)]
            if x[1][-1]==y[1][-1]:
                blocks=defaultdict(Counter)
                for t in q[5:]:blocks[t[:-1]][t[-1]]+=1
                assert all(c[0]==c[210] and c[105]==c[315] for c in blocks.values())
            elif cover(q,420):success=(i,j)
        assert success is not None;witness[h]=[cs[i] for i in success]
    assert len(witness)==6
    return {'profiles':summary,'exceptional_equations':len(exceptions),'uniform_five_row_witnesses':witness}

def literal(row,eq):
    a,b,g,e=row
    for k,l,_,r in eq:
        if k:t=a//k;a-=t*k;b-=t*l;e+=t*r
        else:t=b//l;b-=t*l;e+=t*r
    return (g,a,b,e%840)

def check_refine(directory):
    targets=[list(map(int,l.split())) for l in (directory/'target.txt').read_text().splitlines()]
    assert len(targets)==3
    old=[[(int(g==1),int(g==2),0,e*21) for e,g in zip(t[1:12],t[12:])] for t in targets]
    u=[(0,0,0,21*e) for e in U];one=[(0,0,0,0)]*11
    assert hashlib.sha256((directory/'primary.txt').read_bytes()).digest()==hashlib.sha256((directory/'audit.txt').read_bytes()).digest()
    counts=Counter();n=0
    for l in (directory/'primary.txt').open():
        r=tuple(map(int,l.split()));h=r[0];assert len(r)==68 and 0<=h<3
        eq=[r[j:j+4] for j in (1,5,9) if any(r[j:j+4])]
        x=[(r[j+1],r[j+2],r[j],r[j+4]) for j in range(13,68,5)]
        assert all(r[j+3]==0 and r[j] in (0,1) for j in range(13,68,5))
        counts[len(eq)]+=1;n+=1
        for y in (one,u,old[h]):
            q=[literal(tuple(a-b for a,b in zip(v,w)),eq) for v,w in zip(x,y)]
            assert cover(q,840)
    return n,dict(sorted(counts.items()))

@lru_cache(None)
def solve_pair(e,f):
    a,b,u=e;c,d,v=f;delta=a*d-b*c
    if not delta:return frozenset()
    n=abs(delta);answers=set()
    for i,j in product(range(n),repeat=2):
        x=d*(u+i*L)-b*(v+j*L);y=a*(v+j*L)-c*(u+i*L)
        assert x%delta==0 and y%delta==0
        pt=(x//delta%L,y//delta%L)
        assert (a*pt[0]+b*pt[1]-u)%L==(c*pt[0]+d*pt[1]-v)%L==0
        answers.add(pt)
    assert len(answers)==n
    return frozenset(answers)

def geometry(directory):
    spaces=defaultdict(list);lines=defaultdict(lambda:defaultdict(list));points=defaultdict(lambda:defaultdict(list))
    for i,l in enumerate((directory/'primary.txt').open()):
        a=tuple(map(int,l.split()[:13]));h=a[0]
        eq=[(a[k],a[k+1],4*a[k+3]) for k in (1,5,9) if any(a[k:k+4])]
        if not eq:spaces[h].append(i)
        elif len(eq)==1:
            x,y,r=eq[0];m=gcd(x,y);assert r%m==L%m==0
            for j in range(m):lines[h][x//m,y//m,(r+j*L)//m%L].append(i)
        else:
            for pt in solve_pair(*eq):points[h][pt].append(i)
    expected={}
    for h in spaces:
        planes=spaces[h];curves=lines[h];dots=points[h]
        expected[h,2,0,0,1,0,0,1]=planes
        for (a,b,r),rows in curves.items():
            if a:assert a==1;k=(h,1,r,0,-b,1,0,0)
            else:assert b==1;k=(h,1,0,r,1,0,0,0)
            expected[k]=planes+rows
        for e,f in combinations(curves,2):
            for pt in solve_pair(e,f):dots.setdefault(pt,[])
        for pt,ids in dots.items():
            full=planes+ids
            for (a,b,r),rows in curves.items():
                if (a*pt[0]+b*pt[1]-r)%L==0:full.extend(rows)
            expected[h,0,*pt,0,0,0,0]=full
    sets=set()
    for l in (directory/'strata.txt').open():
        r=tuple(map(int,l.split()));k=r[:8];assert k not in sets and len(r[9:])==r[8]
        assert sorted(expected[k])==sorted(r[9:])
        sets.add(k)
    assert sets==set(expected)
    peaks={}
    for mode in ('primary','audit'):
        graphs=[tuple(map(int,l.split())) for l in (directory/(mode+'_graphs.txt')).read_text().splitlines()]
        assert all(len(r)==12 for r in graphs) and len(graphs)==len(sets) and {r[:8] for r in graphs}==sets
        assert all(r[8]==len(expected[r[:8]]) for r in graphs)
        peaks[mode]=[(r[:8],r[-1]) for r in graphs if r[-1]>6]
    return len(sets),peaks

def reduced(rows):
    def norm(row,mask):
        anchors={};v=[]
        for j,(g,a,e) in enumerate(row):
            g=mask[j]
            if g not in anchors:anchors[g]=a,e
            b,f=anchors[g];v.append((g,a-b,(e-f)%L))
        return tuple(v)
    circles=defaultdict(set)
    for row in rows:
        mask=tuple(x[0] for x in row)
        if any(mask):circles[mask].add(norm(row,mask))
    keep=[c for cs in circles.values() for c in cs]
    for row in rows:
        if any(x[0] for x in row):continue
        if not any(norm(row,mask) in cs for mask,cs in circles.items()):keep.append(row)
    return keep

def hard_lines(directory,failed):
    from circle_bounds import orthogonality_polynomial,polynomial_gcd
    needed=set(k for mode in failed.values() for k,_ in mode)
    assert needed=={(h,1,0,b,1,0,0,0) for h,b in enumerate((1344,2016,2688))}
    ids=set();headers={}
    for l in (directory/'strata.txt').open():
        r=tuple(map(int,l.split()));k=r[:8]
        if k in needed:ids.update(r[9:]);headers[k]=r[9:]
    rows={}
    for i,l in enumerate((directory/'primary.txt').open()):
        if i in ids:
            a=tuple(map(int,l.split()));rows[i]=[a[j:j+5] for j in range(13,68,5)]
    summary={}
    for (h,d,x,y,dx,dy,_,_),indices in headers.items():
        raw={tuple((g,dx*a+dy*b,(4*e+a*x+b*y)%L) for g,a,b,c,e in rows[i]) for i in indices}
        verts=reduced(raw);by=defaultdict(list)
        for v in verts:
            mask=tuple(z[0] for z in v);assert all(e%168==0 for g,a,e in v)
            by[mask].append((mask,tuple(z[1] for z in v),tuple(z[2]*3//2 for z in v)))
        assert sorted(map(len,by.values()))==[6,6,6] and all(any(mask) for mask in by)
        triples=0
        for fs in by.values():
            polys={(i,j):orthogonality_polynomial(fs[i],fs[j]) for i,j in combinations(range(6),2)}
            for i,j,k in combinations(range(6),3):
                divisor=polynomial_gcd(polynomial_gcd(polys[i,j],polys[i,k]),polys[j,k])
                assert len(divisor)==1 and any(divisor[0])
                triples+=1
        summary[h]={'groups':3,'triples_excluded':triples,'maximum_added_row_bound':6}
    return summary

def verify_all():
    if sys.flags.optimize:raise RuntimeError('Run without -O.')
    assert hashlib.sha256((ROOT/'g21_pair_families.txt').read_bytes()).hexdigest()==PAIR_HASH
    for p in (3361,20161):assert all(p%k for k in range(2,int(p**.5)+1))
    rows=census(3250);normals=sorted({normalize(row) for row in rows})
    targets=[list(map(int,l.split())) for l in (ROOT/'3250/target.txt').read_text().splitlines()]
    assert targets==[[i,*k[1::2],*k[::2]] for i,k in enumerate(normals)]
    matches=normal3(census(3297))
    slice_dir=ROOT/'3297/slices'
    before=(slice_dir/'target.txt').read_bytes(),(slice_dir/'mapping.json').read_bytes()
    prepare()
    assert before==((slice_dir/'target.txt').read_bytes(),(slice_dir/'mapping.json').read_bytes())
    mapping=json.loads((ROOT/'3297/slices/mapping.json').read_text());assert mapping['parameter_exponents_mod40']==list(range(2,40,4)) and mapping['normalized_count']==3
    assert len(mapping['assignments'])==60
    result={'normal_forms_3297':{str(k):v for k,v in matches.items()},'generic_3297':verify_generic()}
    for n in ('3250','3297/slices'):
        directory=ROOT/n;count,ranks=check_refine(directory);domains,failed=geometry(directory)
        if n=='3250':assert domains==10251 and not any(failed.values())
        else:assert domains==4335;result['hard_lines']=hard_lines(directory,failed)
        result[n]={'records':count,'ranks':ranks,'domains':domains,'uniform_row_bound':9}
    result['passed']=True
    (ROOT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Both pilots verified, including all phase domains and witnesses.',flush=True)
if __name__=='__main__':verify_all()
