"""Exact integer contradictions for every triple of pair-compatibility loci."""
from pathlib import Path
from collections import defaultdict,Counter
from fractions import Fraction as Q
from itertools import product
from math import gcd,lcm
import json,sys
from verify_rows import cover
R=Path(__file__).resolve().parent;D=R/(sys.argv[1]if len(sys.argv)>1 else'analysis')

def contradiction(rows):
 n=len(rows);basis={}
 for i,row in enumerate(rows):
  a=list(map(Q,row[:3]));w=[Q(j==i)for j in range(n)]
  for pivot,(v,b)in sorted(basis.items()):
   coefficient=a[pivot]
   a=[x-coefficient*y for x,y in zip(a,v)];w=[x-coefficient*y for x,y in zip(w,b)]
  pivot=next((j for j in range(3)if a[j]),None)
  if pivot is None:
   multiple=lcm(*(x.denominator for x in w));weights=[int(multiple*x)for x in w]
   common=gcd(*weights);weights=[x//common for x in weights]
   assert all(sum(v*r[k]for v,r in zip(weights,rows))==0 for k in range(3))
   residue=sum(v*r[3]for v,r in zip(weights,rows))%840
   if residue:return weights,residue
  else:
   coefficient=a[pivot];basis[pivot]=([x/coefficient for x in a],[x/coefficient for x in w])
 return None

def main():
 families=defaultdict(list)
 for l in(D/'point_families.txt').open():
  z=list(map(int,l.split()));families[z[0]].append((z[2:13],z[13:24]))
 pairs=defaultdict(list);triangles={};expected=set()
 for l in(D/'points_graphs.txt').open():
  z=l.split();p,n,m,k=int(z[0]),int(z[2]),int(z[3]),int(z[4]);edges=list(map(int,z[5+k:]));assert len(edges)==2*m
  adj=[set()for _ in range(n)]
  for a,b in zip(edges[::2],edges[1::2]):adj[a].add(b);adj[b].add(a)
  for a in range(n):
   for b in adj[a]:
    if b>a:
     for c in adj[a]&adj[b]:
      if c>b:expected.add((p,a,b,c))
 for l in(D/'triangle_loci.txt').open():
  z=l.split();assert z[0]in('PAIR','TRIANGLE')
  if z[0]=='PAIR':
   key=tuple(map(int,z[1:4]));h=list(map(int,z[4:]));assert len(h)==12
   pairs[key].append(tuple(h))
  else:
   key=tuple(map(int,z[1:5]));stats=tuple(map(int,z[5:]));assert stats[-1]==0 and key not in triangles;triangles[key]=stats
 assert set(triangles)==expected
 pair_records=0
 for (p,a,b),ls in pairs.items():
  assert len(set(ls))==len(ls)
  ga,ea=families[p][a];gb,eb=families[p][b]
  values=[(ga[j],-gb[j],0,(ea[j]-eb[j])%840)for j in range(11)]
  for h in ls:
   constraints=[h[i:i+4]for i in(0,4,8)]
   assert cover(values,constraints),(p,a,b,h)
   pair_records+=1
 checks=0;types=Counter();points=Counter();residues=Counter()
 final=D/'triangle_contradictions.txt';temporary=D/'triangle_contradictions.txt.writing'
 with temporary.open('w')as out:
  for (p,a,b,c),stats in sorted(triangles.items()):
   ab,ac,bc=[sorted(pairs[p,i,j])for i,j in((a,b),(a,c),(b,c))]
   assert stats[:3]==(len(ab),len(ac),len(bc))and stats[3]==len(ab)*len(ac)*len(bc)
   types[sum(max(families[p][i][0])for i in(a,b,c))]+=1;points[p]+=1
   for i,j,k in product(range(len(ab)),range(len(ac)),range(len(bc))):
    equations=[]
    for h,s,t in((ab[i],0,1),(ac[j],0,2),(bc[k],1,2)):
     for start in(0,4,8):
      aa,bb,cc,e=h[start:start+4];assert cc==0
      if not(aa or bb):assert e==0;continue
      row=[0,0,0,e];row[s]=aa;row[t]=bb;equations.append(row)
    cert=contradiction(equations);assert cert is not None,(p,a,b,c,i,j,k,equations)
    w,e=cert;checks+=1;residues[e]+=1
    out.write(' '.join(map(str,[p,a,b,c,i,j,k,len(w)]+w+[e]))+'\n')
 if final.exists():
  assert final.read_bytes()==temporary.read_bytes(),'saved integer contradictions disagree with independent reconstruction'
  temporary.unlink()
 else:temporary.replace(final)
 report={'verified':True,'points_with_triangles':dict(sorted(points.items())),'triangles':len(triangles),
         'triangle_types_by_free_row_circles':dict(sorted(types.items())),'pair_loci_symbolically_verified':pair_records,
         'contradictions':checks,'contradiction_residues':dict(sorted(residues.items())),'feasible_triangles':0}
 (D/'triangle_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
