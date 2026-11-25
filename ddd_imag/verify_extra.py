"""Independent integer checks of symmetries, witnesses and the example."""
from pathlib import Path
from collections import Counter
from itertools import combinations
import json
R=Path(__file__).resolve().parent

def cover(values,m):
 if not values:return True
 anchor=values[0]
 for p in (2,3,5,7):
  if len(values)<p:continue
  rem=list(values)
  for k in range(p):
   v=(anchor+m//p*k)%m
   if v not in rem:break
   rem.remove(v)
  else:
   if cover(rem,m):return True
 return False

def cyclotomic20(exponents):
 # Reduction modulo Phi_20=X^8-X^6+X^4-X^2+1.
 p=[0]*20
 for e in exponents:p[e%20]+=1
 for k in range(19,7,-1):
  v=p[k];p[k]=0
  for d,c in ((0,-1),(2,1),(4,-1),(6,1)):p[k-8+d]+=c*v
 return p[:8]

def example():
 x=[0,12,4,16,3,8,8,13,13,18,3]
 y=[0,16,12,8,0,15,4,19,10,9,5]
 z=[0,16,12,8,0,10,14,4,5,15,4]
 gy=[0,0,0,0,1,1,0,0,1,0,1]
 gz=[0,0,0,0,1,1,0,0,1,1,0]
 original=[list(map(int,l.split()))for l in(R/'reference_R.txt').open()]
 old=[[0]*11]+[[(e//2+a+7*b+2*c)%20 for e,a,b,c in original[i:i+11]]for i in(0,11)]
 for row,groups in ((x,[0]*11),(y,gy),(z,gz)):
  for o in old:
   for g in set(groups):assert cover([(u-v)%20 for u,v,h in zip(row,o,groups)if h==g],20)
 for row,groups in ((y,gy),(z,gz)):
  A=cyclotomic20([u-v for u,v,g in zip(x,row,groups)if not g])
  assert any(A)
  # <x,row(s)>=A+s^{-1}B; B=-eta^19 A.
  assert cyclotomic20([u-v for u,v,g in zip(x,row,groups)if g]+[u-v+19 for u,v,g in zip(x,row,groups)if not g])==[0]*8
  assert cover([(u-v-19*g)%20 for u,v,g in zip(x,row,groups)],20)
 yy=[(e+19*g)%20 for e,g in zip(y,gy)]
 zz=[(e+19*g)%20 for e,g in zip(z,gz)]
 assert cyclotomic20([u-v for u,v in zip(yy,zz)])==[5]+[0]*7
 # Two exact simultaneous y/z choices, showing that this graph edge exists.
 for s,t in ((17,2),(37,22)):
  assert cover([(3*u+s*g-3*v-t*h)%60 for u,g,v,h in zip(y,gy,z,gz)],60)
 data={'point_exponents_mod20':[1,7,2],'fixed_row':x,
       'y_coefficients_mod20':y,'y_groups':gy,'z_coefficients_mod20':z,'z_groups':gz,
       'forced_y_phase_mod20':19,'forced_z_phase_mod20':19,
       'yz_compatible_phase_pairs_mod60':[[17,2],[37,22]],
       'yz_inner_product_at_forced_phases':5,'verified':True}
 (R/'incompatible_triangle_example.json').write_text(json.dumps(data,indent=2)+'\n')
 return data

def main():
 example()
 inp=[list(map(int,l.split()))for l in(R/'b_equals_ic.txt').open()]
 original=[list(map(int,l.split()))for l in(R/'reference_R.txt').open()]
 assert inp==[[(e-10*c)%40,a,b+c,0]for e,a,b,c in original]
 # Sign symmetry on the original three-parameter triple.
 cp=list(range(7))+[8,7,9,10]
 for row in (original[:11],original[11:]):
  for j in range(11):
   e,a,b,c=row[j];assert [(e+20*b)%40,a,b,c]==row[cp[j]]
 # The order-four action (a,b)->(ia,-ib) on this particular surface.
 cp=list(range(7))+[9,10,8,7]
 for row in (inp[:11],inp[11:]):
  for j in range(11):
   e,a,b,c=row[j];assert [(e+10*a-10*b)%40,a,b,c]==row[cp[j]]
 data={}
 for l in(R/'analysis/point_families.txt').open():
  z=list(map(int,l.split()));data.setdefault(z[0],[]).append((z[2:13],z[13:24]))
 seen=set()
 for l in(R/'analysis/points_graphs_audit.txt.witnesses').open():
  i,j,k,s,t=map(int,l.split());assert i not in seen;seen.add(i);a,b=divmod(i,120)
  rows=[[0]*11]+[[(42*v[0]+14*a*v[1]+14*b*v[2])%1680 for v in inp[q:q+11]]for q in(0,11)]
  for index,phase in ((j,s),(k,t)):
   groups,e=data[i][index];rows.append([(2*x+g*phase)%1680 for g,x in zip(groups,e)])
  for u,v in combinations(rows,2):assert cover([(x-y)%1680 for x,y in zip(u,v)],1680),(i,j,k)
 assert seen==set(range(14400))
 result={'verified':True,'input_specialization':True,'sign_symmetry':True,
         'order_four_point_symmetry':True,'incompatible_triangle_example':True,
         'point_fifth_row_witnesses':14400}
 (R/'extra_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
if __name__=='__main__':main()
