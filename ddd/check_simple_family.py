"""Verify the three-parameter displayed family symbolically, with integer phases."""
from pathlib import Path
import json

u=[0,4,8,12,16,0,10,5,15,0,10]
v=[0,8,16,4,12,0,5,10,15,10,0]
ub=[[0,0,0] for _ in range(5)]+[[1,0,0] for _ in range(4)]+[[1,0,1] for _ in range(2)]
vb=[[0,0,0] for _ in range(5)]+[[0,1,0] for _ in range(4)]+[[0,1,1] for _ in range(2)]
w=[(a-b)%20 for a,b in zip(u,v)];wb=[[a-b for a,b in zip(x,y)] for x,y in zip(ub,vb)]
supports=[[[0,1,2,3,4],[5,6],[7,8],[9,10]],[[0,1,2,3,4],[5,7],[6,8],[9,10]],[[0,1,2,3,4],[5,9],[6,7],[8,10]]]
for c,B,P in zip((u,v,w),(ub,vb,wb),supports):
    for block in P:
        j=block[0];p=len(block)
        assert sorted((c[k]-c[j])%20 for k in block)==[20*k//p for k in range(p)]
        assert all(B[k]==B[j] for k in block)
root=Path(__file__).resolve().parent
(root/'simple_family.json').write_text(json.dumps({'modulus':20,'exponents':[u,v],'basis':[ub,vb],'supports':supports,'verified':True},indent=2)+'\n')
print('Displayed three-parameter family verified symbolically.')
