from collections import Counter
from pathlib import Path
from functools import lru_cache
import sys,json
M=840

def reduce(v,rows):
 v=list(v)
 for i,r in enumerate(rows):
  if not r[i]:continue
  q=v[i]//r[i]
  for k in range(i,3):v[k]-=q*r[k]
  v[3]+=q*r[3]
 v[3]%=M
 return tuple(v)

def minus(a,b):return tuple(x-y for x,y in zip(a,b))
def cover(values,rows):
 if not values:return True
 a=values[0];left=[reduce(minus(x,a),rows)for x in values]
 for p in (2,3,5,7):
  if len(left)<p:continue
  rem=left[:]
  for k in range(p):
   v=(0,0,0,M//p*k)
   if v not in rem:break
   rem.remove(v)
  else:
   if cover(rem,rows):return True
 return False

def readrow(path):
 z=[list(map(int,l.split()))for l in Path(path).read_text().splitlines()]
 return [[(0,0,0,0)]*11]+[[tuple(v[1:]+[21*v[0]])for v in z[i:i+11]]for i in(0,11)]
def unpack(z):
 constraints=[z[i:i+4]for i in(0,4,8)];groups=[z[12+5*j]for j in range(11)];values=[tuple(z[13+5*j:17+5*j])for j in range(11)]
 return constraints,groups,values

def verify(file,oldfile):
 old=readrow(oldfile);n=0;counts=Counter();generic=[]
 for line in Path(file).open():
  z=list(map(int,line.split()));assert len(z)==67;r,g,v=unpack(z)
  for o in old:
   for group in set(g):assert cover([minus(v[j],o[j])for j in range(11)if g[j]==group],r),(n,group)
  rank=sum(bool(r[i][i])for i in range(3));counts[(rank,max(g))]+=1
  if rank==0:generic.append({'groups':g,'values':v})
  n+=1
 print('symbolically verified',n,'rows; (parameter constraints,additional row parameters):',dict(counts),flush=True)
 print('generic',len(generic),flush=True)
 edges=[]
 for i in range(len(generic)):
  for j in range(i):
   a,b=generic[i],generic[j]
   if max(a['groups'])or max(b['groups']):continue
   if cover([minus(x,y)for x,y in zip(a['values'],b['values'])],[[0]*4]*3):edges.append((j,i))
 print('uniform generic compatible pairs',edges,flush=True)
 return {'rows':n,'counts':{str(k):v for k,v in counts.items()},'generic':generic,'generic_edges':edges}
if __name__=='__main__':
 ans=verify(sys.argv[1],sys.argv[2])
 if len(sys.argv)>3:Path(sys.argv[3]).write_text(json.dumps(ans,indent=2)+'\n')
