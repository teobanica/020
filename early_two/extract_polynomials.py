"""Exact bivariate phase matrices and minor polynomials over Q(zeta_60)."""
from pathlib import Path
from itertools import permutations,combinations
import json
from cyclotomic60 import K,POW
ROOT=Path(__file__).resolve().parent
SEED=json.loads((ROOT/'seed.json').read_text())['exponents']

def det(A):
 n=len(A);s=K()
 for p in permutations(range(n)):
  v=K((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n)))
  for i,j in enumerate(p):v=v*A[i][j]
  s=s+v
 return s

def build(line,H):
 w=line.split();b=list(map(int,w[1:12]));e=list(map(int,w[12:]));k=max(b)+1;assert all(v%7==0 for v in e)
 matrices=[]
 for mask in range(4):
  A=[[K() for _ in range(k)] for r in range(4)]
  for r in range(4):
   for j in range(11):
    if r==2 and bool(mask&1)!=(j>=5):continue
    if r==3 and bool(mask&2)!=(j>=5):continue
    A[r][b[j]]=A[r][b[j]]+POW[(e[j]//7-3*H[r+1][j])%60]
  matrices.append(A)
 return k,matrices

def polynomials(k,matrices):
 out=[]
 if k<=4:
  for rows in combinations(range(4),k):
   dq=int(2 in rows);ds=int(3 in rows);P={}
   for i in range(dq+1):
    for j in range(ds+1):
     A=matrices[i+2*j];v=det([[A[r][c] for c in range(k)] for r in rows])
     if v:P[i,j]=v
   out.append(P)
 for rows in combinations(range(4),k-1):
  dq=int(2 in rows);ds=int(3 in rows);N=[]
  for omit in range(k):
   co={}
   for i in range(dq+1):
    for j in range(ds+1):
     A=matrices[i+2*j];co[i,j]=det([[A[r][c] for c in range(k) if c!=omit] for r in rows])
   P={}
   for (i,j),a in co.items():
    for (u,v),b in co.items():P[i-u+dq,j-v+ds]=P.get((i-u+dq,j-v+ds),K())+a*b.conj()
   N.append(P)
  for P in N[1:]:out.append({e:P[e]-N[0][e] for e in P if P[e]!=N[0][e]})
 return out
