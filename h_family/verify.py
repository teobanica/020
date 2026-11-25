"""Independent integer checks of seed, census and every saved point rejection.
The full driver is required to reproduce the large exhaustive search itself.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations,product
from collections import Counter
import json,math
ROOT=Path(__file__).resolve().parent
L=352800
COUNTS={'B':712800,'C':443520,'D':166320,'E':123200,'F':34650}
SIZES={'B':(7,2,2),'C':(5,3,3),'D':(5,2,2,2),'E':(3,3,3,2),'F':(3,2,2,2,2)}
H=json.loads((ROOT/'seed.json').read_text())['exponents']
def stat(text):
    words=text.split();out={};i=0
    while i<len(words):
        if words[i]=='dimensions':out['dimensions']=list(map(int,words[i+1:]));break
        out[words[i]]=int(words[i+1]);i+=2
    return out
def regular(x,r):
    y=[(x[j]-H[r][j]*(L//20))%L for j in range(11)]
    common=[];specific={}
    for prime in (2,3,5,7):
        wanted=[h*(L//prime)for h in range(prime)]
        for support in combinations(range(11),prime):
            left=[j for j in support if j<5];right=[j for j in support if j>=5]
            mask=sum(1<<j for j in support)
            if not left or not right:
                if sorted((y[j]-y[support[0]])%L for j in support)==wanted:common.append(mask)
                continue
            j,k=left[0],right[0]
            for h in range(1,prime):
                q=(y[j]-y[k]-h*(L//prime))%L
                if sorted((y[col]+(q if col>=5 else 0)-y[j])%L for col in support)==wanted:specific.setdefault(q,[]).append(mask)
    for extra in [[],*specific.values()]:
        cycles=common+extra;bycolumn=[[mask for mask in cycles if mask>>j&1]for j in range(11)]
        @lru_cache(None)
        def cover(mask):
            if not mask:return True
            anchor=(mask&-mask).bit_length()-1
            return any(cover(mask^cycle) for cycle in bycolumn[anchor] if mask&cycle==cycle)
        if cover(2047):return True
    return False
def main():
    seed=json.loads((ROOT/'seed.json').read_text());assert len(H)==5 and H[0]==[0]*11
    assert {tuple(x['rows'])for x in seed['quotient_certificates']}==set(combinations(range(5),2))
    for rec in seed['quotient_certificates']:
        a,b=rec['rows'];v=[(H[a][j]-H[b][j])%20 for j in range(11)]
        blocks=[rec['pentagon']]+rec['pairs'];assert sorted(j for s in blocks for j in s)==list(range(11))
        assert all(sorted((v[j]-v[s[0]])%20 for j in s)==[h*(20//len(s))for h in range(len(s))]for s in blocks)
        assert all(all(j<5 for j in s)or all(j>=5 for j in s)for s in blocks)
    for c,sizes in SIZES.items():
        expected=math.factorial(11)//(math.prod(sizes)*math.prod(math.factorial(n)for n in Counter(sizes).values()))
        assert expected==COUNTS[c]
    report={'presentations':sum(COUNTS.values()),'seed_quotients':10,'arithmetic':{},'point_rejections':{}}
    cached={}
    for mode in ('primary','audit'):
        records=json.loads((ROOT/(mode+'_summary.json')).read_text());assert {v['profile']for v in records}==set(COUNTS)
        report['arithmetic'][mode]={}
        for record in records:
            c=record['profile'];d=ROOT/mode/c
            stats={name:stat(record[name])for name in ('search','refine','bivariate')}
            first=stats['search'];ref=stats['refine'];last=stats['bivariate']
            assert first['tests']==COUNTS[c] and sum(first['dimensions'])==first['leaves']
            assert first['free']==ref['free']==first['actual']==ref['actual']==last['actual']==0
            assert first['families']==ref['records']==len((d/'families.txt').read_text().splitlines())
            assert ref['unresolved']==last['count']==len((d/'refine_unresolved.txt').read_text().splitlines())
            assert last['unresolved']==0 and not(d/'bivariate_unresolved.txt').read_text()
            point_count=0
            for filename,key,n in [('points.txt','search','passed'),('refine_points.txt','refine','passed'),('bivariate_points.txt','bivariate','normpassed')]:
                lines=(d/filename).read_text().splitlines();assert len(lines)==stats[key][n]
                for line in lines:
                    v=list(map(int,line.split()));assert len(v)==18 and v[13]==0
                    x=tuple(v[2:13]);assert all(0<=e<L for e in x)
                    if x not in cached:cached[x]=next((r for r in range(1,5)if not regular(x,r)),None)
                    assert cached[x] is not None
                point_count+=len(lines)
            report['point_rejections'][mode+'_'+c]=point_count
            report['arithmetic'][mode][c]=stats
    for c in COUNTS:
        a=report['arithmetic']['primary'][c]['search'];b=report['arithmetic']['audit'][c]['search']
        for key in ('tests','leaves','dimensions'):assert a[key]==b[key],(c,key)
    report['distinct_point_rejections']=len(cached)
    report['unresolved']=0;report['regular_sixth_rows']=0
    (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items()if k!='arithmetic'},indent=2))
if __name__=='__main__':main()
