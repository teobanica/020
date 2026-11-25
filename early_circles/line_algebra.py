from pathlib import Path
from itertools import permutations,combinations
import json,time
from cyclotomic60 import K,POW,rref
ROOT=Path(__file__).resolve().parent
def det(A):
 n=len(A);out=K()
 for s in permutations(range(n)):
  v=K((-1)**sum(s[i]>s[j] for i in range(n) for j in range(i+1,n)))
  for i,j in enumerate(s):v=v*A[i][j]
  out=out+v
 return out
def trim(p):
 p=list(p)
 while p and not p[-1]:p.pop()
 return p
def divmodp(a,b):
 a=trim(a);b=trim(b);q=[K()] * max(0,len(a)-len(b)+1)
 if not b:raise ValueError()
 inv=b[-1].inverse()
 while len(a)>=len(b):
  d=len(a)-len(b);v=a[-1]*inv;q[d]=v
  for j,w in enumerate(b):a[j+d]=a[j+d]-v*w
  a=trim(a)
 return trim(q),a
def gcd(a,b):
 while b:a,b=b,divmodp(a,b)[1]
 if not a:return []
 inv=a[-1].inverse();return [v*inv for v in a]
def make(line,H):
 w=line.split();b=list(map(int,w[1:12]));e=list(map(int,w[12:]));k=max(b)+1
 assert all(v%7==0 for v in e)
 A0=[[K() for _ in range(k)] for _ in range(4)];A1=[[K() for _ in range(k)] for _ in range(4)]
 for r in range(4):
  for j in range(11):
   mat=A1 if r==3 and j>=5 else A0
   mat[r][b[j]]=mat[r][b[j]]+POW[(e[j]//7-3*H[r+1][j])%60]
 for r in range(3):A1[r]=A0[r][:]
 return k,A0,A1
def minor(A,rs,cs):return det([[A[r][c] for c in cs] for r in rs])
def polynomials(k,A0,A1):
 if k<=4:
  for rs in combinations(range(4),k):
   a=minor(A0,rs,range(k));b=minor(A1,rs,range(k)) if 3 in rs else K()
   yield trim([a,b])
 for rs in combinations(range(4),k-1):
  pol=[]
  for omit in range(k):
   cs=[j for j in range(k) if j!=omit];a=minor(A0,rs,cs);b=minor(A1,rs,cs) if 3 in rs else K()
   pol.append([a*b.conj(),a*a.conj()+b*b.conj(),b*a.conj()])
  for p in pol[1:]:yield trim([v-w for v,w in zip(p,pol[0])])
