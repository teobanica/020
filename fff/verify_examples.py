"""Exact symbolic checks of the two displayed examples (no numerical tolerance)."""
from pathlib import Path
from math import gcd
import json
ROOT=Path(__file__).resolve().parent

def main():
    # Monomials: a sixth-root exponent and six integer parameter exponents.
    zero=(0,)*6
    def m(c,k=None):
        a=[0]*6
        if k is not None:a[k]=1
        return c%6,tuple(a)
    # Parameters (a,b,s,t,p,q).
    H=[[m(0) for _ in range(11)] for _ in range(3)]
    H[1][:3]=[m(0),m(2),m(4)]
    H[2][:3]=[m(0),m(4),m(2)]
    for start,aa,ss,tt in [(3,0,2,3),(7,1,4,5)]:
        for j,(c,uses_a) in enumerate([(0,False),(3,False),(0,True),(3,True)]):
            powers=[0]*6;powers[ss]=1;powers[aa]=int(uses_a)
            H[1][start+j]=(c,tuple(powers))
        for j,c in enumerate([0,0,3,3]):
            powers=[0]*6;powers[tt]=1
            H[2][start+j]=(c,tuple(powers))
    certificates=[]
    for i in range(3):
        for j in range(i+1,3):
            row=[((a[0]-b[0])%6,tuple(x-y for x,y in zip(a[1],b[1])))
                 for a,b in zip(H[i],H[j])]
            assert sorted(x[0] for x in row[:3])==[0,2,4]
            assert all(x[1]==zero for x in row[:3])
            pairs=[]
            for start in [3,7]:
                left=list(range(start,start+4))
                while left:
                    a=left.pop(0)
                    b=next(x for x in left if (row[x][0]-row[a][0])%6==3
                           and row[x][1]==row[a][1])
                    left.remove(b);pairs.append([a,b])
            certificates.append({"rows":[i,j],"triangle":[0,1,2],"pairs":pairs})
    # Six parameter exponent vectors span Q^6: recover s,t,p,q directly,
    # then a,b by ratios within the second row.
    recovered=[H[1][3][1],H[2][3][1],H[1][7][1],H[2][7][1],
               tuple(x-y for x,y in zip(H[1][5][1],H[1][3][1])),
               tuple(x-y for x,y in zip(H[1][9][1],H[1][7][1]))]
    assert set(recovered)=={tuple(int(i==j) for i in range(6)) for j in range(6)}
    ex=json.loads((ROOT/"fff_summary.json").read_text())["examples"]["r24"]
    U,V=ex["exponents"]
    assert U==[0,16,8,10,22,18,6,15,3,6,18]
    assert V==[0,16,22,8,10,0,12,3,0,15,12]
    for row,part in zip([U,V,[(a-b)%24 for a,b in zip(U,V)]],ex["supports"]):
        for block in part:
            assert sorted((row[k]-row[block[0]])%24 for k in block)==list(range(0,24,24//len(block)))
    assert gcd(24,U[3]-U[7])==1
    expected=[[int(j in [5,6]),int(j in [7,8,9,10])] for j in range(11)]
    assert ex["basis"]==[expected,expected]
    report={"verified":True,"six_parameter_family":{"modulus":6,"monomials":H,
            "quotient_certificates":certificates,"dimension":6},
            "order24_example":ex,"primitive_cross_ratio_exponent":-5,
            "fourth_row_obstruction_scope":"only retaining the common F3 block"}
    (ROOT/"example_verification.json").write_text(json.dumps(report,indent=2)+"\n")
    print("Verified the six-parameter family and the order-24 example.")

if __name__=="__main__":
    main()

