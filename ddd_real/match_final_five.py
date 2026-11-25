from pathlib import Path
from itertools import permutations
import json
R=Path(__file__).resolve().parent;M=840
j=json.loads((R/'analysis/generic_published.json').read_text())
z=[list(map(int,l.split()))for l in(R/'b_equals_c.txt').read_text().splitlines()]
old=[[(0,0,0,0)]*11]+[[tuple(v[1:]+[21*v[0]])for v in z[i:i+11]]for i in(0,11)]
E=[[0]*11,[0,12,24,36,48,3,33,30,0,0,30],[0,24,48,12,51,36,21,0,30,0,30],[0,36,12,48,24,24,24,44,44,4,4],[0,48,36,24,12,12,12,52,52,32,32]]
T=[[14*x for x in row]for row in E]

def minus(a,b):return tuple(a[k]-b[k]for k in range(3))+((a[3]-b[3])%M,)
def match(F):
 for perm in permutations(range(5)):
  A=[[minus(F[perm[r]][col],F[perm[0]][col])for col in range(11)]for r in range(5)]
  # Desired zero-character rows 0,3,4; singleton rows 1,2 have one tail character.
  if any(any(v[:3]!=(0,0,0)for v in A[r])for r in (0,3,4)):continue
  if any(any(v[:3]!=(0,0,0)for v in A[r][:7])or len({v[:3]for v in A[r][7:]})!=1 for r in (1,2)):continue
  u=A[1][7][:2];v=A[2][7][:2]
  if abs(u[0]*v[1]-u[1]*v[0])!=1:continue
  for anchor in range(7):
   B=[[minus(v,row[anchor])for v in row]for row in A]
   lookup={tuple(B[r][col][3]for r in range(5)):col for col in range(7)}
   target=[tuple(T[r][col]for r in range(5))for col in range(7)]
   if any(k not in lookup for k in target):continue
   first=[lookup[k]for k in target]
   if len(set(first))!=7:continue
   for last in permutations(range(7,11)):
    if any(B[r][col][3]!=T[r][7+k]for r in (0,3,4)for k,col in enumerate(last)):continue
    d=[[ (B[r][col][3]-T[r][7+k])%M for k,col in enumerate(last)]for r in (1,2)]
    if any(row[0]!=row[1]or row[2]!=row[3]for row in d):continue
    delta=(d[0][2]-d[0][0])%M
    if delta!=(d[1][2]-d[1][0])%M:continue
    return {'row_permutation':perm,'column_permutation':first+list(last),'anchor':anchor,'A_powers':u,'B_powers':v,'A_root':d[0][2],'B_root':d[1][0],'C_over_B_root':delta}
 return None
out=[]
for edge in j['generic_edges']:
 F=old+[[tuple(v)for v in j['generic'][i]['values']]for i in edge]
 witness=match(F);out.append({'edge':edge,'witness':witness});print(edge,witness,flush=True)
(R/'analysis/five_equivalences.json').write_text(json.dumps(out,indent=2)+'\n')
