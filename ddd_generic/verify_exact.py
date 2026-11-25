"""Independent integer checks for the generic DDD triple calculation."""
from pathlib import Path
from itertools import combinations
from collections import Counter
from functools import lru_cache
import json,hashlib,math
ROOT=Path(__file__).resolve().parent
M=840

def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def cycles(values,columns=None):
 if not values:return []
 if columns is None:columns=list(range(len(values)))
 a=values[0];diff=[sub(x,a)for x in values];diff=[v[:3]+(v[3]%M,)for v in diff]
 for p in (2,3,5,7):
  if p>len(values):continue
  ids=[]
  for k in range(p):
   value=(0,0,0,M//p*k)
   try:i=next(i for i,v in enumerate(diff)if i not in ids and v==value)
   except StopIteration:break
   ids.append(i)
  else:
   rest=[i for i in range(len(values))if i not in ids]
   rec=cycles([values[i]for i in rest],[columns[i]for i in rest])
   if rec is not None:return [[columns[i]for i in ids]]+rec
 return None

def polynomial_division(a,b):
 a=a[:]
 while a and not a[-1]:a.pop()
 q=[0]*max(0,len(a)-len(b)+1)
 while len(a)>=len(b):
  assert b[-1]==1
  d=len(a)-len(b);c=a[-1];q[d]=c
  for j,x in enumerate(b):a[d+j]-=c*x
  while a and not a[-1]:a.pop()
 return q,a
@lru_cache(None)
def cyclotomic(n):
 a=[-1]+[0]*(n-1)+[1]
 for d in range(1,n):
  if n%d==0:
   a,r=polynomial_division(a,cyclotomic(d));assert not r
 return a

def main():
 inp=ROOT/'family_3286.txt'
 if not inp.exists():inp=ROOT.parent/'work_ddd_positive/family_3286.txt'
 z=[list(map(int,l.split()))for l in inp.read_text().splitlines()];assert len(z)==22
 old=[[(0,0,0,0)]*11]+[[tuple(x[1:]+[21*x[0]])for x in z[i:i+11]]for i in(0,11)]
 for r in old:
  for k in range(3):assert max(v[k]for v in r)-min(v[k]for v in r)<=1
 for a,b in combinations(old,2):assert cycles([sub(x,y)for x,y in zip(a,b)])is not None
 allrows=set();profiles={};stats={}
 for p,total in zip('BCDEF',(712800,443520,166320,123200,34650)):
  primary=ROOT/f'3286_{p}.txt';audit=ROOT/f'audit_3286_{p}.txt'
  s=set(primary.read_text().splitlines());t=set(audit.read_text().splitlines());assert s==t,(p,len(s),len(t));allrows|=s;profiles[p]=s
  for tag in ('3286','audit_3286'):
   logs=(ROOT/f'{tag}_{p}.log').read_text().splitlines();summary=next(l for l in logs if l.startswith('tests '));words=summary.split();assert int(words[1])==total
   assert int(words[5])==len(s)
   assert int(words[3])==dict(B=0,C=276,D=1800,E=0,F=11264)[p]
  stats[p]={'presentations':total,'fourth_rows':len(s),'sha256':hashlib.sha256(primary.read_bytes()).hexdigest()}
 rows=[]
 for l in sorted(allrows):
  k=list(map(int,l.split()));assert len(k)==67 and k[:12]==[0]*12
  g=k[12::5];v=[tuple(k[13+5*j:17+5*j])for j in range(11)];assert g==[0]*11
  assert all(vj[:3]==(0,0,0)for vj in v)
  certs=[cycles([sub(x,y)for x,y in zip(v,r)])for r in old];assert all(c is not None for c in certs)
  assert all(x[3]%7==0 for x in v)
  rows.append({'exponents_mod120':[x[3]//7 for x in v],'cycles_against_old':certs,'profile':'C'if l in profiles['C']else'D'})
 assert len(rows)==8 and len(profiles['C'])==len(profiles['D'])==4
 edges=[];orthogonal=[];pair_certificates={};phi=cyclotomic(120)
 for a,b in combinations(range(8),2):
  e=[(x-y)%120 for x,y in zip(rows[a]['exponents_mod120'],rows[b]['exponents_mod120'])]
  poly=[0]*120
  for x in e:poly[x]+=1
  if not polynomial_division(poly,phi)[1]:orthogonal.append([a,b])
  cert=cycles([(0,0,0,7*x)for x in e])
  if cert is not None:edges.append([a,b]);pair_certificates[f'{a},{b}']=cert
 assert edges==[[0,4],[1,5]]
 assert orthogonal==edges,(orthogonal,edges)
 # Swap the two last column pairs and apply (a,b,c)->(-ab/c,c,b).
 perm=list(range(7))+[9,10,7,8];T=((1,1,-1),(0,0,1),(0,1,0));shift=(M//2,0,0)
 for r in old:
  for j in range(11):
   v=r[j];new=tuple(sum(v[k]*T[k][d]for k in range(3))for d in range(3))+( (v[3]+sum(v[k]*shift[k]for k in range(3)))%M,)
   assert new==r[perm[j]],(j,new,r[perm[j]])
 for a,b in((0,1),(4,5)):
  assert [rows[a]['exponents_mod120'][j]for j in perm]==rows[b]['exponents_mod120']
 # The existing G rows s,t are exactly rows 1,5 of this list.
 assert rows[1]['exponents_mod120']==[0,72,24,96,48,48,48,88,88,8,8]
 assert rows[5]['exponents_mod120']==[0,96,72,48,24,24,24,104,104,64,64]
 exceptional=set()
 for p in 'BCDEF':
  for l in (ROOT/f'audit_3286_{p}.txt.exceptional').read_text().splitlines():
   a,b,c,e=map(int,l.split());assert 0<max(abs(a),abs(b),abs(c))<=5 and 0<=e<840 and e%2==0;exceptional.add((a,b,c,e))
 bychar={}
 for a,b,c,e in exceptional:
  d=math.gcd(math.gcd(abs(a),abs(b)),abs(c));v=(a//d,b//d,c//d)
  bychar[v]=math.lcm(bychar.get(v,1),420*d)
 expected={(0,0,1):420,(0,1,-1):840,(0,1,0):420,(1,-1,0):420,(1,-1,1):420,(1,0,-1):420,(1,0,0):840,(1,1,-2):420,(1,1,-1):840,(1,1,0):420,(1,2,-2):420,(2,1,-2):420,(2,1,-1):840,(3,1,-1):420,(3,2,-2):420}
 assert bychar==expected and len(exceptional)==5001
 ans={'verified':True,'first_presentations':sum(s['presentations']for s in stats.values()),'profiles':stats,'rows':rows,'edges':edges,'orthogonal_pairs':orthogonal,'pair_certificates':pair_certificates,'exceptional_equations':len(exceptional),'character_cover':[{'powers':list(v),'root_order':m}for v,m in sorted(bychar.items())]}
 (ROOT/'verification.json').write_text(json.dumps(ans,indent=2)+'\n')
 (ROOT/'exceptional_equations.txt').write_text(''.join(' '.join(map(str,v))+'\n'for v in sorted(exceptional)))
 print(json.dumps({k:ans[k]for k in ('verified','first_presentations','edges','orthogonal_pairs','exceptional_equations','character_cover')},indent=2))
if __name__=='__main__':main()
