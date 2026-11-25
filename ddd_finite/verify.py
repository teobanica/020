"""Independent integer checks of all retained affine rows and graph certificates."""
from pathlib import Path
from functools import lru_cache
from collections import Counter,defaultdict
from itertools import combinations
import json,math,hashlib,time
ROOT=Path(__file__).resolve().parent

def determinant(matrix):
 a=[r[:]for r in matrix];n=len(a);den=1;sign=1
 for k in range(n-1):
  if not a[k][k]:
   q=next((j for j in range(k+1,n)if a[j][k]),None)
   if q is None:return 0
   a[k],a[q]=a[q],a[k];sign=-sign
  pivot=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    num=pivot*a[i][j]-a[i][k]*a[k][j];assert num%den==0;a[i][j]=num//den
  den=pivot
  for i in range(k+1,n):a[i][k]=0
 return sign*a[-1][-1]
def multiply(a,b):return [[sum(x*y for x,y in zip(row,col))for col in zip(*b)]for row in a]
@lru_cache(maxsize=400000)
def cycles(values,modulus):
 if not values:return True
 anchor=values[0]
 for p in (2,3,5,7):
  if p>len(values)or modulus%p:continue
  left=list(values)
  for h in range(p):
   x=(anchor+modulus//p*h)%modulus
   if x not in left:break
   left.remove(x)
  else:
   if cycles(tuple(left),modulus):return True
 return False
def regular(values,modulus):
 values=sorted(values)
 if not values:return True
 # All common root rotations describe the same cycle-cover problem.
 key=min(tuple(sorted((x-a)%modulus for x in values))for a in set(values))
 return cycles(key,modulus)
def symbolic(family,old,modulus=840):
 c,groups,e=family
 for group in range(c):
  values=[(e[j]-(modulus//40)*old[j])%modulus for j in range(11)if groups[j]==group]
  if not regular(values,modulus):return False
 return True
def stat(line):
 z=line.split();return {z[i]:int(z[i+1])for i in range(0,z.index('dimensions')if 'dimensions'in z else len(z),2)}
def main():
 start=time.monotonic();records=json.loads((ROOT/'classes.json').read_text());assert len(records)==452;classes={r['class_id']:r for r in records}
 if(ROOT/'audit_label_hashes.json').exists():
  for name,digest in json.loads((ROOT/'audit_label_hashes.json').read_text()).items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
 assert Counter(r['dimension']for r in records)=={0:218,1:190,2:41,3:3}
 for r in records:
  A,U,V=r['A'],r['U'],r['V'];assert abs(determinant(U))==abs(determinant(V))==1
  D=multiply(multiply(U,A),V);assert all(D[i][j]==(r['diagonal'][i]if i==j and i<r['rank']else 0)for i in range(7)for j in range(6))
  assert r['dimension']==6-r['rank'];assert all(sum(U[i][j]*r['g'][j]for j in range(7))%2==0 for i in range(r['rank'],7))
 count=0
 for line in(ROOT/'families.txt').open():
  z=list(map(int,line.split()));idx,d,*_=z[:6];r=classes[idx];assert d==r['dimension'] and 1<=d<=3 and len(z)==28+22*d
  e=[[0]*11,z[6:17],z[17:28]];a=[[(0,)*d for j in range(11)],[],[]];pos=28
  for row in range(1,3):
   for j in range(11):a[row].append(tuple(z[pos:pos+d]));pos+=d
  for (left,right),support in zip(((1,0),(2,0),(1,2)),r['supports']):
   for block in support:
    anchor=block[0];power=tuple(a[left][anchor][k]-a[right][anchor][k]for k in range(d));base=e[left][anchor]-e[right][anchor]
    assert all(tuple(a[left][j][k]-a[right][j][k]for k in range(d))==power for j in block)
    assert sorted((e[left][j]-e[right][j]-base)%40 for j in block)==[40//len(block)*k for k in range(len(block))]
  count+=1
 assert count==201248;print('positive-dimensional affine identities',count,flush=True)
 original=[]
 for line in(ROOT/'zero_seeds.txt').open():
  z=list(map(int,line.split()));assert len(z)==23+z[22];key=tuple(z[:22]);original.append(key)
  assert all(0<=v<40 for v in key)
  u=key[::2];v=key[1::2]
  for q in (u,v,tuple((x-y)%40 for x,y in zip(u,v))):assert regular(q,40)
 assert len(original)==len(set(original))==54836
 indexed={};exception=set()
 for line in(ROOT/'zero_index.txt').open():
  z=list(map(int,line.split()));idx,n=z[:2];supports=z[2:2+n];uv=z[2+n:];assert len(uv)==22 and idx not in indexed
  indexed[idx]=[[0]*11,uv[:11],uv[11:]]
  if 826 in supports:exception.add(idx)
 assert set(indexed)==set(range(54836)) and len(exception)==632
 example=json.loads((ROOT/'new_five_example.json').read_text());ex=example['exponents'];assert example['modulus']==120 and example['seed']in exception
 assert ex[:3]==[[3*x for x in row]for row in indexed[example['seed']]]
 pairs=set()
 for item in example['cycles']:
  a,b=item['rows'];pairs.add((a,b));assert sorted(j for block in item['cycles']for j in block)==list(range(11))
  for block in item['cycles']:
   p=len(block);assert p in(2,3,5);base=ex[a][block[0]]-ex[b][block[0]]
   assert [(ex[a][j]-ex[b][j]-base)%120 for j in block]==[120//p*k for k in range(p)]
 assert pairs==set(combinations(range(5),2)) and math.gcd(ex[3][2],120)==1
 labeltotals=defaultdict(lambda:[0,0,0]);labelids=set()
 for line in(ROOT/'labels_counts.txt').read_text().splitlines():
  idx,dim,branches,labels,components,finite=map(int,line.split());r=classes[idx];assert idx not in labelids and dim==r['dimension'] and branches==math.prod(r['diagonal']) and components==labels*branches;labelids.add(idx)
  labeltotals[dim][0]+=1;labeltotals[dim][1]+=labels;labeltotals[dim][2]+=components
 assert labelids==set(classes) and dict(labeltotals)=={0:[218,576960,1227520],1:[190,94768,194336],2:[41,3168,6480],3:[3,216,432]}

 f=(ROOT/'zero_groups.txt').open();ng=int(next(f));assert ng==36;group_seeds={};allids=[]
 for _ in range(ng):
  header=list(map(int,next(f).split()));g,n=header[:2];u=header[2:];ids=[]
  for j in range(n):
   z=list(map(int,next(f).split()));idx=z[0];assert indexed[idx][1:]==[u,z[1:]];ids.append(idx)
  group_seeds[g]=ids;allids.extend(ids)
 assert sorted(allids)==list(range(54836))
 report={'support_classes':452,'zero_dimensional_support_classes':218,'finite_triples':54836,'positive_affine_components_verified':count,'groups':36,'heights':Counter(),'exception_heights':Counter(),'fourth_row_records':Counter(),'seeds_with_circles':0,'edge_witnesses':0,'edges':0,'max_graph_clique':0,'per_group':[]}
 # The packaged audit stores byte-identical duplicate files by their hashes.
 hashes=json.loads((ROOT/'audit_hashes.json').read_text())if(ROOT/'audit_hashes.json').exists()else None
 for g in range(36):
  folder=ROOT/'zero_results'/str(g);audit=ROOT/'zero_audit'/str(g)
  for name in ('pair_families.txt','pair_summary.txt','extensions.txt'):
   digest=hashlib.sha256((folder/name).read_bytes()).hexdigest()
   if(audit/name).exists():assert digest==hashlib.sha256((audit/name).read_bytes()).hexdigest(),(g,name)
   else:assert hashes is not None and digest==hashes[str(g)][name]
  pair=(folder/'pair_summary.txt').read_text().splitlines();assert len(pair)==5
  expected=[712800,443520,166320,123200,34650]
  total=0
  for line,n in zip(pair,expected):
   z=line.split();assert int(z[z.index('tests')+1])==n;total+=int(z[z.index('unique')+1])
  assert sum(1 for line in(folder/'pair_families.txt').open())==total
  oldrows=group_seeds[g];data=defaultdict(list)
  for line in(folder/'extensions.txt').open():
   z=list(map(int,line.split()));assert len(z)==24;idx,c=z[:2];groups,e=z[2:13],z[13:];assert idx in oldrows and c in(1,2) and set(groups)==set(range(c))
   assert all(0<=x<840 for x in e)
   for k in range(c):assert e[groups.index(k)]==0
   family=(c,groups,e)
   assert all(symbolic(family,old)for old in indexed[idx]),(g,idx,family)
   data[idx].append(family);report['fourth_row_records'][c-1]+=1
  summary=(folder/'extension_summary.txt').read_text().splitlines()[-1].split();assert summary[:2]==['summary','group']
  assert int(summary[summary.index('seeds')+1])==len(oldrows)
  assert int(summary[summary.index('extension_families')+1])==sum(map(len,data.values()))
  assert int(summary[summary.index('extendible')+1])==len(data)
  for idx,rows in data.items():
   assert len(rows)==len({(c,tuple(gs),tuple(e))for c,gs,e in rows})
   report['seeds_with_circles']+=any(c==2 for c,gs,e in rows)
  graphfile=folder/'possible_graphs.txt';auditgraph=folder/'possible_graphs_audit.txt';assert graphfile.read_bytes()==auditgraph.read_bytes()
  witnesses={}
  for line in Path(str(auditgraph)+'.witnesses').read_text().splitlines():
   idx,a,b,q,s=map(int,line.split());assert idx not in witnesses;witnesses[idx]=(a,b,q,s)
  seen=set();edgegraphs=set();gh=Counter()
  for line in graphfile.read_text().splitlines():
   z=line.split();idx=int(z[0]);kind=z[1];n,m,omega=map(int,z[2:5]);best=list(map(int,z[5:5+omega]));edges=list(map(int,z[5+omega:]));assert len(edges)==2*m and n==len(data[idx]) and idx not in seen;seen.add(idx)
   assert kind==('PARAMETRIC'if any(c==2 for c,gs,e in data[idx])else'FINITE')
   adjacency=[set()for j in range(n)];edge_set=set()
   for a,b in zip(edges[::2],edges[1::2]):assert 0<=a<b<n and(a,b)not in edge_set;edge_set.add((a,b));adjacency[a].add(b);adjacency[b].add(a)
   for a,b in edge_set:assert not(adjacency[a]&adjacency[b]),(idx,a,b)
   assert omega==(2 if m else 1) and len(best)==omega and len(set(best))==omega and all(0<=x<n for x in best)
   if omega==2:assert tuple(sorted(best))in edge_set
   if m:
    edgegraphs.add(idx);a,b,q,s=witnesses[idx];assert(a,b)in edge_set
    rows=[[42*x for x in old]for old in indexed[idx]]
    for k,t in ((a,q),(b,s)):
     c,gs,e=data[idx][k];rows.append([(2*e[j]+gs[j]*t)%1680 for j in range(11)])
    for a,b in combinations(range(5),2):assert regular([(x-y)%1680 for x,y in zip(rows[a],rows[b])],1680),(idx,a,b)
    report['edge_witnesses']+=1
   height=3+omega;report['heights'][height]+=1;gh[height]+=1
   if idx in exception:report['exception_heights'][height]+=1
   report['edges']+=m;report['max_graph_clique']=max(report['max_graph_clique'],omega)
  assert seen==set(data) and edgegraphs==set(witnesses)
  for idx in set(oldrows)-seen:
   report['heights'][3]+=1;gh[3]+=1
   if idx in exception:report['exception_heights'][3]+=1
  report['per_group'].append({'group':g,'seed_count':len(oldrows),'pair_families':total,'heights':dict(gh)})
  print('verified group',g,dict(gh),flush=True)
 assert sum(report['heights'].values())==54836 and sum(report['exception_heights'].values())==632 and report['max_graph_clique']==2
 report['seconds']=round(time.monotonic()-start,3);report['triangle_free']=True
 (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items()if k!='per_group'},indent=2))
if __name__=='__main__':main()
