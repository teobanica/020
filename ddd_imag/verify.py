"""Verify the complete retained imaginary-surface certificate exactly."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import gcd
import hashlib,json
from verify_rows import verify as symbolic_verify,readrow,cover,minus
from verify_extra import main as extra_verify
from verify_triangles import main as triangle_verify
R=Path(__file__).resolve().parent

def token_rows(path):return {tuple(map(int,l.split()))for l in path.open()}
def graphs(path,audit,curve=False):
 assert path.read_bytes()==audit.read_bytes()
 seen=set();hist=Counter();cliques=Counter();triangle_points={};maxn=0
 for l in path.open():
  z=l.split()
  if curve:
   i,t,n,m,ta,tb,tc=map(int,z[:7]);key=(i,t);assert(ta,tb,tc)==(-1,-1,-1);ed=list(map(int,z[7:]));assert t==-1 or 0<=t<840 and t%7
  else:
   i=int(z[0]);key=i;assert z[1]in('FINITE','PARAMETRIC');n,m,k=map(int,z[2:5]);assert k in(2,3);cliques[k]+=1;ed=list(map(int,z[5+k:]))
  assert key not in seen;seen.add(key);assert len(ed)==2*m
  adj=[set()for _ in range(n)]
  for a,b in zip(ed[::2],ed[1::2]):assert 0<=a<b<n;adj[a].add(b);adj[b].add(a)
  assert sum(map(len,adj))==2*m
  triangles={(a,b,c)for a,b in zip(ed[::2],ed[1::2])for c in adj[a]&adj[b]if c>b}
  if curve:assert not triangles
  else:
   assert bool(triangles)==(k==3)
   if triangles:triangle_points[i]=len(triangles)
   for a,b,c in triangles:assert not(adj[a]&adj[b]&adj[c])
  hist[m]+=1;maxn=max(maxn,n)
 if curve:
  nc=json.loads((R/'analysis/preparation.json').read_text())['curve_count']
  assert seen=={(i,t)for i in range(nc)for t in range(-1,840)if t==-1 or t%7}
 else:assert seen==set(range(14400))
 return {'graphs':len(seen),'maximum_vertices':maxn,'edge_histogram':dict(sorted(hist.items())),
         'clique_histogram':dict(sorted(cliques.items())),'points_with_triangles':triangle_points,
         'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def main():
 D=R/'analysis';allrows=set();profiles={}
 old=readrow(R/'b_equals_ic.txt')
 for u,v in combinations(old,2):assert cover([minus(x,y)for x,y in zip(u,v)],[[0]*4]*3)
 for p,total in zip('BCDEF',(712800,443520,166320,123200,34650)):
  a=token_rows(R/f'primary_{p}.txt');b=token_rows(R/f'audit_{p}.txt');assert a==b;allrows|=a
  sums=[]
  for tag in('primary','audit'):
   line=next(l for l in(R/f'{tag}_{p}.log').read_text().splitlines()if l.startswith('tests '));v=line.split()
   assert int(v[1])==total and int(v[5])==len(a);sums.append(v)
  assert sums[0]==sums[1]
  profiles[p]={'first_presentations':total,'intermediate_records':int(v[3]),'final_records':len(a)}
 assert allrows==token_rows(R/'all_families.txt')
 points=set();curves=set();hist=Counter();diagonals=Counter()
 for r in allrows:
  assert len(r)==67 and r[10]==r[2]==r[6]==0
  rank=sum(bool(r[i])for i in(0,5,10));hist[rank,max(r[12::5])]+=1
  if rank==0:assert max(r[12::5])==0 and all(r[j]==0 for j in range(13,67)if (j-13)%5<3)
  elif rank==1:
   a,b,c,e=r[:4]if r[0]else r[4:8];d=gcd(abs(a),abs(b));assert max(r[12::5])==0
   for k in range(d):
    q=(Q(e,840)+k)/d%1;assert 60%q.denominator==0;curves.add((a//d,b//d,q))
  else:
   a,b,c,e=r[:4];u,v,w,f=r[4:8];det=a*v-b*u;assert det>0;pts=set();diagonals[a,v]+=1
   for k in range(det):
    for ell in range(det):
     aa=Q(v*(e+840*k)-b*(f+840*ell),840*det)%1;bb=Q(-u*(e+840*k)+a*(f+840*ell),840*det)%1
     pts.add((aa,bb));assert 120%aa.denominator==0 and 120%bb.denominator==0
   assert len(pts)==det;points|=pts;assert max(r[12::5])<=1
 intersections=set()
 for (a,b,q),(c,d,s)in combinations(curves,2):
  det=a*d-b*c
  if not det:continue
  pts=set()
  for k in range(abs(det)):
   for ell in range(abs(det)):
    aa=(d*(q+k)-b*(s+ell))/det%1;bb=(-c*(q+k)+a*(s+ell))/det%1
    pts.add((aa,bb));assert 120%aa.denominator==0 and 120%bb.denominator==0
  assert len(pts)==abs(det);intersections|=pts
 prep=json.loads((D/'preparation.json').read_text());assert not prep['partial']
 assert len(curves)==prep['curve_count'] and len(points)==prep['rank_two_support_points'] and len(intersections)==prep['curve_intersections']
 assert {(v['normal'][0],v['normal'][1],Q(v['root_mod840'],840))for v in prep['curves']}==curves
 # The complementary generic graph is constant in both old parameters.
 sym=symbolic_verify(R/'all_families.txt',R/'b_equals_ic.txt')
 assert len(sym['generic'])==8 and len(sym['generic_edges'])==2
 assert sym['rows']==len(allrows) and sym['counts']=={str(k):v for k,v in hist.items()}
 adj=[set()for _ in range(8)]
 for a,b in sym['generic_edges']:adj[a].add(b);adj[b].add(a)
 for a,b in sym['generic_edges']:assert not(adj[a]&adj[b])
 (D/'generic_verification.json').write_text(json.dumps(sym,indent=2)+'\n')
 points_report=graphs(D/'points_graphs.txt',D/'points_graphs_audit.txt')
 curves_report=graphs(D/'curves_graphs.txt',D/'curves_graphs_audit.txt',True)
 assert (D/'triangle_loci.txt').read_bytes()==(D/'triangle_loci_audit.txt').read_bytes()
 triangle_verify()
 tri=json.loads((D/'triangle_verification.json').read_text());assert tri['feasible_triangles']==0
 assert {int(k):v for k,v in tri['points_with_triangles'].items()}==points_report['points_with_triangles']
 extra_verify()
 report={'verified':True,'profiles':profiles,'records':len(allrows),
  'family_histogram':{str(k):v for k,v in sorted(hist.items())},
  'rank_two_diagonals':{str(k):v for k,v in sorted(diagonals.items())},
  'rank_one_connected_curves':len(curves),'rank_two_points':len(points),
  'curve_intersections':len(intersections),'point_graphs':points_report,'curve_graphs':curves_report,
  'generic_rows':8,'generic_edges':2,'simultaneous_phase_check':tri,'regular_six_row_extensions':0}
 (R/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items()if k!='simultaneous_phase_check'},indent=2))
if __name__=='__main__':main()
