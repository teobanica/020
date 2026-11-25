from pathlib import Path
import json
from cyclotomic60 import K,POW,rref
def solve():
    seed=json.loads(Path('seed.json').read_text());H=seed['exponents'];out=[]
    for idx,line in enumerate(Path('unresolved.txt').read_text().splitlines()):
        w=line.split();b=list(map(int,w[1:12]));e=list(map(int,w[12:]));assert all(x%7==0 for x in e)
        A=[[sum((POW[(e[j]//7-3*H[r][j])%60] for j in range(11) if b[j]==k),K()) for k in range(5)] for r in range(1,5)]
        R,piv=rref(A);free=[j for j in range(5) if j not in piv];assert len(free)==2
        N=[[K(j==f) for f in free] for j in range(5)]
        for i,j in enumerate(piv):N[j]=[-R[i][f] for f in free]
        assert all(sum((A[r][j]*N[j][t] for j in range(5)),K())==0 for r in range(4) for t in range(2))
        rec={'index':idx,'rank':len(piv),'free_columns':free,'kernel_basis':[[x.serial() for x in row] for row in N]}
        zero=next((j for j,row in enumerate(N) if not any(row)),None)
        if zero is not None:
            rec['exclusion']='forced_zero_cycle_phase';rec['column']=zero
        else:
            # Free coordinates s,t must both be unit. Divide by s and put z=t/s.
            # |a+bz|^2=1 becomes conj(a)b*z+a*conj(b)*conj(z)=1-|a|^2-|b|^2.
            eq=[[a.conj()*b,a*b.conj(),1-a*a.conj()-b*b.conj()] for a,b in N]
            augmented=[row+[K(i==j) for j in range(5)] for i,row in enumerate(eq)]
            rr,pp=rref(augmented,2)
            assert all(sum((row[3+i]*eq[i][j] for i in range(5)),K())==row[j] for row in rr for j in range(3))
            inconsistent=any(not any(row[:2]) and row[2] for row in rr)
            if inconsistent:
                rec['exclusion']='inconsistent_linear_norm_equations'
                rec['norm_equations']=[[x.serial() for x in row] for row in eq]
                witness=next(row for row in rr if not any(row[:2]) and row[2])
                rec['inconsistency_combination']=[x.serial() for x in witness[3:]]
                rec['nonzero_constant']=witness[2].serial()
            elif len(pp)==2:
                u,v=rr[0][2],rr[1][2]
                assert u*v!=1 or v!=u.conj()
                rec['exclusion']='unique_norm_candidate_not_unit'
                rec['z']=u.serial();rec['zbar']=v.serial();rec['product_minus_one']=(u*v-1).serial()
            else:
                raise AssertionError(('unresolved exact norm system',idx,len(pp)))
        out.append(rec);print(idx,rec['exclusion'],flush=True)
    Path('exact_remaining.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':solve()
