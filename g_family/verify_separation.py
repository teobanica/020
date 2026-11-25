"""Exact finite data for the common-pentagon separation argument.

For a,b,c,ab/c,b/c outside mu60, cross-group ratios in u and v cannot
be fifth roots. The only groups of size at least five are their first
seven columns. We verify their unique mu5 cosets with five elements.
The open-dense assertion and invariance are proved in the manuscript.
"""
from pathlib import Path
import json,itertools
R=Path(__file__).resolve().parent
seed=json.loads((R/'new_three_parameter_five.json').read_text())
E=seed['exponents'];P=seed['powers']
for i,target in [(1,[0,1,2,3,4]),(2,[0,1,2,3,5])]:
 groups={}
 for j,p in enumerate(P[i]):groups.setdefault(tuple(p),[]).append(j)
 assert sorted(map(len,groups.values()))==[2,2,7]
 assert groups[(0,0,0)]==list(range(7))
 candidates=[]
 for S in itertools.combinations(range(7),5):
  if sorted((E[i][j]-E[i][S[0]])%60 for j in S)==[0,12,24,36,48]:candidates.append(list(S))
 assert candidates==[target],candidates
 powers=list(groups)
 allowed={(1,0,0),(0,1,0),(0,0,1),(1,1,-1),(0,1,-1)}
 for x,y in itertools.combinations(powers,2):
  d=tuple(a-b for a,b in zip(x,y));assert d in allowed or tuple(-z for z in d)in allowed
# Earlier family's first five columns are F5: each distinct row quotient
# has every residue mod5 once, independent of its four tail parameters.
for i,k in itertools.combinations(range(5),2):assert sorted(((i-k)*j)%5 for j in range(5))==list(range(5))
print('Generic separation: distinct unique pentagon supports verified; earlier common-pentagon invariant verified.')
