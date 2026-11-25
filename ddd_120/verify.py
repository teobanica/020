"""Independent exact checks of the finite-level DDD certificate."""
from pathlib import Path
from collections import Counter,defaultdict
from itertools import combinations
from functools import lru_cache
from math import gcd
import json,hashlib
from verify_rows import verify as symbolic_verify,readrow,cover as symbolic_cover,minus,reduce
R=Path(__file__).resolve().parent;M=120

def token_rows(path):return {tuple(map(int,l.split()))for l in path.open()}

@lru_cache(maxsize=None)
def cycle_cover(values):
 if not values:return ()
 anchor=values[0]
 for p in(5,3,2):
  if p>len(values):continue
  rem=list(values);cycle=tuple((anchor+M//p*k)%M for k in range(p))
  for x in cycle:
   if x not in rem:break
   rem.remove(x)
  else:
   found=cycle_cover(tuple(rem))
   if found is not None:return (cycle,)+found
 return None

@lru_cache(maxsize=None)
def regular(values):
 # A quick necessary homomorphism, then an exact multiset cycle partition.
 # ROOT has order 120 modulo 241; this identity is checked at entry.
 if sum(POWERS[e]for e in values)%241:return False
 return cycle_cover(values)is not None

ROOT=next(z for z in range(2,241)if pow(z,M,241)==1 and all(pow(z,M//p,241)!=1 for p in(2,3,5)))
POWERS=tuple(pow(ROOT,e,241)for e in range(M))
def quotient(a,b):return tuple(sorted((x-y)%M for x,y in zip(a,b)))

def clique_exists(adjacency,candidates,need):
 if need==0:return True
 while candidates.bit_count()>=need:
  bit=candidates&-candidates;v=bit.bit_length()-1;candidates^=bit
  if clique_exists(adjacency,candidates&adjacency[v],need-1):return True
 return False

def main():
 D=R/'analysis';profiles={};allrows=set()
 original=[list(map(int,l.split()))for l in(R/'reference_R.txt').open()]
 inp=[list(map(int,l.split()))for l in(R/'ratio_zeta3.txt').open()]
 assert inp==[[(3*e-40*c)%M,a,b+c,0]for e,a,b,c in original]
 # R(a,-b,c) equals R(a,b,c) with columns 7 and 8 exchanged.
 cp=list(range(7))+[8,7,9,10]
 for row in(original[:11],original[11:]):
  for j in range(11):
   e,a,b,c=row[j];assert [(e+20*b)%40,a,b,c]==row[cp[j]]
 # Swapping the final column pairs sends (a,b,c) to (-ab/c,c,b).
 cp=list(range(7))+[9,10,7,8]
 for row in(original[:11],original[11:]):
  for j in range(11):
   e,a,b,c=row[j];assert [(e+20*a)%40,a,a+c,b-a]==row[cp[j]]
 oldformal=readrow(R/'ratio_zeta3.txt')
 for u,v in combinations(oldformal,2):assert symbolic_cover([minus(x,y)for x,y in zip(u,v)],[[0]*4]*3)
 generic=[[0,72,24,96,48,48,48,8,8,88,88],[0,96,72,48,24,24,24,64,64,104,104]]
 for row in generic:
  for o in oldformal:assert symbolic_cover([minus((0,0,0,e),x)for e,x in zip(row,o)],[[0]*4]*3)
 assert regular(quotient(*generic))
 for p,total in zip('CDEF',(443520,166320,123200,34650)):
  a=token_rows(R/f'primary_{p}.txt');b=token_rows(R/f'audit_{p}.txt');assert a==b;allrows|=a
  counts=[]
  for tag in('primary','audit'):
   line=next(l for l in(R/f'{tag}_{p}.log').read_text().splitlines()if l.startswith('tests '));v=line.split()
   assert int(v[1])==total and int(v[5])==len(a);counts.append(tuple(map(int,v[1:6:2])))
  assert counts[0]==counts[1]
  profiles[p]={'first_presentations':total,'intermediate_records':counts[0][1],'final_records':len(a)}
 assert allrows==token_rows(D/'all_families.txt')
 hist=Counter();diagonals=Counter()
 for r in allrows:
  constraints=[r[i:i+4]for i in(0,4,8)]
  assert r[0]>0 and r[5]>0 and all(r[i]==0 for i in(2,4,6,8,9,10,11))
  assert 0<=r[1]<r[5] and M%r[0]==M%r[5]==0
  # The complete old-parameter locus lies in mu_120 squared.
  assert reduce((M,0,0,0),constraints)==(0,0,0,0)
  assert reduce((0,M,0,0),constraints)==(0,0,0,0)
  hist[max(r[12::5])]+=1;diagonals[r[0],r[5]]+=1
 assert max(hist)<=1
 sym=symbolic_verify(D/'all_families.txt',R/'ratio_zeta3.txt')
 (D/'symbolic_verification.json').write_text(json.dumps(sym,indent=2)+'\n')
 assert sym['rows']==len(allrows)
 assert (D/'rows.txt').read_bytes()==(D/'rows_audit.txt').read_bytes()
 assert (D/'graphs.txt').read_bytes()==(D/'graphs_audit.txt').read_bytes()
 rows=defaultdict(list)
 for l in(D/'rows.txt').open():
  z=tuple(map(int,l.split()));assert len(z)==12 and z[1]==0 and all(0<=x<M for x in z[1:]);rows[z[0]].append(z[1:])
 assert set(rows)==set(range(M*M))
 nhist=Counter();ehist=Counter();khist=Counter();seen=set();examples={};total_pairs=0;level120_example=None
 for l in(D/'graphs.txt').open():
  z=list(map(int,l.split()));point,n,m,k=z[:4];best=z[4:4+k];ev=z[4+k:]
  assert point not in seen;seen.add(point);assert len(ev)==2*m and n==len(rows[point])
  assert rows[point]==sorted(set(rows[point]));adj=[0]*n;edges=set(zip(ev[::2],ev[1::2]));assert len(edges)==m
  for i,j in edges:assert 0<=i<j<n;adj[i]|=1<<j;adj[j]|=1<<i
  assert len(set(best))==k and all(0<=i<n for i in best)
  for i,j in combinations(best,2):assert adj[i]>>j&1
  assert not clique_exists(adj,(1<<n)-1,k+1),(point,k)
  a,b=divmod(point,M);old=[[0]*11]+[[(e+a*x+b*y)%M for e,x,y,_ in inp[q:q+11]]for q in(0,11)]
  assert all(tuple(row)in rows[point]for row in generic)
  for row in rows[point]:
   for o in old:assert regular(quotient(row,o)),(point,row,o)
  for j in range(n):
   for i in range(j):assert regular(quotient(rows[point][i],rows[point][j]))==((i,j)in edges),(point,i,j)
  nhist[n]+=1;ehist[m]+=1;khist[k]+=1;total_pairs+=n*(n-1)//2
  if k not in examples:
   matrix=old+[list(rows[point][i])for i in best]
   decomposition=[]
   for i,j in combinations(range(len(matrix)),2):
    cycles=cycle_cover(quotient(matrix[i],matrix[j]));assert cycles is not None
    decomposition.append({'row_pair':[i,j],'cycles_exponents_mod120':cycles})
   d=M
   for row in matrix:
    for e in row:d=gcd(d,e)
   examples[k]={'a_exponent':a,'b_exponent':b,'c_exponent':(b-40)%M,'matrix_exponents_mod120':matrix,
                'minimal_Butson_level':M//d,'pairwise_cycle_partitions':decomposition}
  if point==120:
   matrix=old+generic;decomposition=[]
   for i,j in combinations(range(5),2):
    cycles=cycle_cover(quotient(matrix[i],matrix[j]));assert cycles is not None
    decomposition.append({'row_pair':[i,j],'cycles_exponents_mod120':cycles})
   assert matrix[1][9]==1
   level120_example={'a_exponent':1,'b_exponent':0,'c_exponent':80,'matrix_exponents_mod120':matrix,
                     'minimal_Butson_level':120,'pairwise_cycle_partitions':decomposition}
 assert seen==set(range(M*M))
 prep=json.loads((D/'preparation.json').read_text());assert not prep['partial'] and prep['records']==len(allrows)
 result={'verified':True,'modulus':M,'profiles':profiles,'records':len(allrows),
   'record_diagonals':{str(k):v for k,v in sorted(diagonals.items())},'free_row_parameters':dict(sorted(hist.items())),
   'point_graphs':len(seen),'actual_rows':sum(n*c for n,c in nhist.items()),'edges':sum(n*c for n,c in ehist.items()),
   'pair_tests':total_pairs,'maximum_vertices':max(nhist),'row_count_histogram':dict(sorted(nhist.items())),
   'edge_histogram':dict(sorted(ehist.items())),'maximum_added_rows_histogram':dict(sorted(khist.items())),
   'maximum_total_rows':3+max(khist),'input_specialization_verified':True,'sign_symmetry_verified':True,'inverse_symmetry_verified':True,
   'point_family_count':prep['point_family_total'],'points_with_free_rows':len(prep['points_with_free_rows']),
   'root_valued_examples':examples,'level120_example':level120_example,'uniform_fifth_row_witness':generic,
   'scope':'The labelled R(a,b,c), b/c in {zeta_3,zeta_3^(-1),-zeta_3,-zeta_3^(-1)}, a,b in mu_120; all added rows in mu_120.'}
 (R/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items()if k not in('root_valued_examples','level120_example','edge_histogram')},indent=2),flush=True)
if __name__=='__main__':main()
