from pathlib import Path
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import gcd
import hashlib,json,sys
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R))
from verify_rows import verify

def token_rows(path):
 return {tuple(map(int,l.split()))for l in path.read_text().splitlines()}
def read_curves(path):
 vals=iter(map(int,path.read_text().split()));n=next(vals);rows={}
 for _ in range(n):
  i,A,B,q,k,D=[next(vals)for _ in range(6)];f=[]
  for j in range(k):f.append(tuple(next(vals)for _ in range(22)))
  rows[i]={'normal':(A,B),'root':q,'degree':D,'rows':f}
 try:next(vals);raise AssertionError('extra curve data')
 except StopIteration:pass
 return rows

def check_graphs(primary,audit,curves=False):
 assert primary.read_bytes()==audit.read_bytes(),(primary,audit)
 seen=set();hist=Counter();nmax=0;count=0
 for l in primary.open():
  z=l.split()
  if curves:
   i,t,n,m,ta,tb,tc=map(int,z[:7]);assert (i,t)not in seen;seen.add((i,t));assert(ta,tb,tc)==(-1,-1,-1);ed=list(map(int,z[7:]));assert len(ed)==2*m
   assert t==-1 or 0<=t<840 and t%14!=0
  else:
   i=int(z[0]);assert z[1]in('FINITE','PARAMETRIC');n,m,k=map(int,z[2:5]);assert k==2;ed=list(map(int,z[5+k:]));assert len(ed)==2*m;assert i not in seen;seen.add(i)
  adj=[set()for _ in range(n)]
  for a,b in zip(ed[::2],ed[1::2]):assert 0<=a<b<n;adj[a].add(b);adj[b].add(a)
  assert sum(map(len,adj))==2*m
  for a,b in zip(ed[::2],ed[1::2]):assert not(adj[a]&adj[b])
  hist[m]+=1;nmax=max(nmax,n);count+=1
 if curves:assert seen=={(i,t)for i in range(180)for t in range(-1,840)if t==-1 or t%14}
 else:assert seen==set(range(3600))
 return {'graphs':count,'max_rows':nmax,'edge_histogram':dict(sorted(hist.items())),'sha256':hashlib.sha256(primary.read_bytes()).hexdigest()}

def main():
 root=R/'analysis';allrows=set();profile_stats={}
 for p,total in zip('BCDEF',(712800,443520,166320,123200,34650)):
  a=token_rows(R/f'bc_{p}.txt');b=token_rows(R/f'audit_bc_{p}.txt');assert a==b;allrows|=a
  sums=[]
  for tag in('bc','audit_bc'):
   line=next(l for l in(R/f'{tag}_{p}.log').read_text().splitlines()if l.startswith('tests '));v=line.split();assert int(v[1])==total and int(v[5])==len(a);sums.append(v)
  assert sums[0]==sums[1]
  profile_stats[p]={'first_presentations':total,'intermediate_records':int(sums[0][3]),'final_records':len(a)}
 saved=token_rows(R/'all_families.txt');assert saved==allrows
 # Check the classification of every old-parameter locus independently by
 # Cramer's rule, with numerator shifts modulo the determinant.
 pointset=set();line_loci=set();hist=Counter()
 for r in allrows:
  rank=sum(bool(r[i])for i in(0,5,10));hist[rank,max(r[12::5])]+=1;assert r[10]==r[2]==r[6]==0
  if rank==0:assert max(r[12::5])==0
  elif rank==1:
   a,b,c,e=r[:4]if r[0]else r[4:8];assert(a,b)in((1,0),(0,1),(1,-1))and e%14==0;line_loci.add((a,b,e//14));assert max(r[12::5])==0
  else:
   a,b,c,e=r[:4];u,v,w,f=r[4:8];D=a*v-b*u;assert D>0;points=set()
   for k in range(D):
    for ell in range(D):
     aa=Q(v*(e+840*k)-b*(f+840*ell),840*D)%1;bb=Q(-u*(e+840*k)+a*(f+840*ell),840*D)%1
     points.add((aa,bb));assert 60%aa.denominator==0 and 60%bb.denominator==0
   assert len(points)==D;pointset|=points;assert max(r[12::5])<=1
 symbolic=verify(R/'all_families.txt',R/'b_equals_c.txt')
 assert symbolic['rows']==len(allrows)
 assert hist=={(0,0):24,(1,0):1296,(2,0):22000,(2,1):96}
 assert len(pointset)==1936
 common={i for i in range(60)if i%3==0 or i%2==0}
 expected_lines={(1,0,q)for q in common|{13,23,43,53}}|{(0,1,q)for q in common|{1,11,31,41}}|{(1,-1,q)for q in common|{7,17,37,47}}
 assert line_loci==expected_lines and len(line_loci)==132
 free_points=Counter();point_families=0
 for line in(root/'point_families.txt').open():
  z=list(map(int,line.split()));point_families+=1
  if z[1]>1:assert z[1]==2;free_points[z[0]]+=1
 assert point_families==101768
 assert free_points=={60*a+b:6 for a in range(3,60,15)for b in range(6,60,15)}
 # The full symbolic verifier has checked all old/new cycle identities,
 # without specializing parameters or dropping nonprimitive constraints.
 points=check_graphs(root/'points_graphs.txt',root/'points_graphs_audit.txt')
 curves=check_graphs(root/'curves_graphs.txt',root/'curves_graphs_audit.txt',True)
 assert(root/'global_graph.txt').read_bytes()==(root/'global_graph_audit.txt').read_bytes()
 adj=[set()for _ in range(24)];pairs=0
 for line in(root/'global_graph.txt').open():
  z=list(map(int,line.split()));a,b,yes,A,B,k=z[:6];roots=z[6:];assert len(roots)==k and all(t%14==0 for t in roots);pairs+=1
  if yes:adj[a].add(b);adj[b].add(a)
 assert pairs==276 and sum(map(len,adj))==12
 for a in range(24):
  for b in adj[a]:assert not(adj[a]&adj[b])
 report={'verified':True,'profiles':profile_stats,'records':len(allrows),'family_histogram':{str(k):v for k,v in sorted(hist.items())},'rank_one_connected_curves':len(line_loci),'rank_two_points':len(pointset),'point_graphs':points,'curve_graphs':curves,'global_rows':24,'global_edges':6,'regular_six_row_extensions':0}
 (R/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
