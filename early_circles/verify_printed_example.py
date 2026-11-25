"""Direct verification of the exact sixth row printed in Section 12.3."""
from pathlib import Path
import json
from cyclotomic60 import POW
from quadratic11 import L
ROOT=Path(__file__).resolve().parent
def verify():
 I=POW[15];d=L(0,1);q=(L(-6+3*I)+d*(1+2*I))/10;a=(L(-7)+d*(5*I))/18
 assert q*q.conj()==1 and a*a.conj()==1
 eta=POW[12];x=[L(1),a*eta,a*eta*eta,a*eta*eta*eta,a*eta*eta*eta*eta,a,-a*I,L(I),a*I,L(-I),L(-1)]
 seed=json.loads((ROOT/'seed.json').read_text())['exponents'];t=[L(1),q.conj(),L(1),L(-I),L(-1)]
 old=[[L(POW[3*seed[r][j]])*(t[r] if j>=5 else 1) for j in range(11)] for r in range(5)]
 assert all(sum((x[j]*row[j].conj() for j in range(11)),L())==0 for row in old)
 covers={0:[[1,2,3,4,5],[0,10],[6,8],[7,9]],2:[[1,2,3,4,5],[0,9],[6,8],[7,10]],3:[[1,2,3,4,8],[0,10],[5,6],[7,9]],4:[[1,2,3,4,6],[0,7],[5,8],[9,10]]}
 for r,blocks in covers.items():
  w=[x[j]*old[r][j].conj() for j in range(11)]
  assert sorted(j for block in blocks for j in block)==list(range(11))
  for block in blocks:
   values=[w[j] for j in block];p=len(block);assert all(values[0]*POW[60*k//p] in values for k in range(p))
 w=[x[j]*old[1][j].conj() for j in range(11)]
 assert w==[L(1),a,a,a,a,a*q,a*q*I,-q,a*q,-q,-q*I]
 values=[L(1),a,a*q,a*q*I,-q,-q*I]
 assert all(v+u for j,v in enumerate(values) for u in values[j:])
 return True
if __name__=='__main__':assert verify();print('Printed sixth row and all displayed identities verified exactly.')
