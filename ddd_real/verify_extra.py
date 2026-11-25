from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from math import gcd
import json
R=Path(__file__).resolve().parent;M=840

def minus(x,y):return(x[0]-y[0],x[1]-y[1],(x[2]-y[2])%M)
def cover(values,m):
 if not values:return True
 anchor=values[0]
 for p in(2,3,5,7):
  if len(values)<p:continue
  rem=list(values)
  for k in range(p):
   v=(anchor+m//p*k)%m
   if v not in rem:break
   rem.remove(v)
  else:
   if cover(rem,m):return True
 return False

def main():
 # A different component decomposition: exact rational angles, no search
 # over residues in a presumed containing root group.
 union=set()
 for l in(R/'exceptional_equations_previous.txt').open():
  a,b,c,e=map(int,l.split());d=gcd(gcd(abs(a),abs(b)),abs(c));v=(a//d,b//d,c//d)
  for k in range(d):
   q=(Q(e,840)+k)/d;q%=1;assert 840%q.denominator==0;union.add((v,int(q*840)))
 report=json.loads((R/'connected_tori.json').read_text());saved={(tuple(t['normal']),e)for t in report['tori']for e in t['roots_mod840']};assert union==saved and len(union)==5016
 # Independently substitute every ordinary-equivalence witness into G.
 inp=[list(map(int,l.split()))for l in(R/'b_equals_c.txt').open()]
 old=[[(0,0,0)]*11]+[[(v[1],v[2],21*v[0])for v in inp[i:i+11]]for i in(0,11)]
 published=json.loads((R/'analysis/generic_published.json').read_text());equiv=json.loads((R/'analysis/five_equivalences.json').read_text());G=json.loads((R/'reference_G.json').read_text());assert len(equiv)==6
 assert {tuple(x['edge'])for x in equiv}=={tuple(x)for x in published['generic_edges']}
 for item in equiv:
  w=item['witness'];assert w is not None;F=old+[[tuple(v[:2]+[v[3]])for v in published['generic'][i]['values']]for i in item['edge']]
  rp=w['row_permutation'];cp=w['column_permutation'];assert sorted(rp)==list(range(5))and sorted(cp)==list(range(11))
  A=tuple(w['A_powers'])+(w['A_root'],);B=tuple(w['B_powers'])+(w['B_root'],);C=B[:2]+((B[2]+w['C_over_B_root'])%M,);params=(A,B,C)
  assert abs(A[0]*B[1]-A[1]*B[0])==1
  for i in range(5):
   for j in range(11):
    source=minus(minus(F[rp[i]][cp[j]],F[rp[0]][cp[j]]),minus(F[rp[i]][cp[0]],F[rp[0]][cp[0]]))
    power=G['powers'][i][j];target=(sum(power[k]*params[k][0]for k in range(3)),sum(power[k]*params[k][1]for k in range(3)),(14*G['exponents'][i][j]+sum(power[k]*params[k][2]for k in range(3)))%M)
    assert source==target,(i,j,source,target)
 # Sign symmetry for the original three-parameter triple.
 orig=[list(map(int,l.split()))for l in(R/'reference_R.txt').open()];perm=list(range(7))+[8,7,9,10]
 for row in (orig[:11],orig[11:]):
  for j in range(11):
   e,a,b,c=row[j];f=row[perm[j]];assert [int((e+20*b)%40),a,b,c]==f
 # Every point has a saved actual five-row witness, including row circles.
 data={}
 for l in(R/'analysis/point_families.txt').open():
  z=list(map(int,l.split()));data.setdefault(z[0],[]).append((z[2:13],z[13:24]))
 seen=set()
 for l in(R/'analysis/points_graphs_audit.txt.witnesses').open():
  i,j,k,s,t=map(int,l.split());assert i not in seen;seen.add(i);a,b=divmod(i,60)
  rows=[[0]*11]+[[(42*v[0]+28*a*v[1]+28*b*v[2])%1680 for v in inp[q:q+11]]for q in(0,11)]
  for index,phase in ((j,s),(k,t)):
   groups,e=data[i][index];rows.append([(2*x+g*phase)%1680 for g,x in zip(groups,e)])
  for u,v in combinations(rows,2):assert cover([(x-y)%1680 for x,y in zip(u,v)],1680),(i,j,k)
 assert seen==set(range(3600))
 result={'connected_tori':5016,'five_family_equivalences':6,'entry_identities':330,'sign_symmetry':True,'point_fifth_row_witnesses':3600,'verified':True}
 (R/'extra_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
if __name__=='__main__':main()
