"""Exact Q(zeta_60) elimination of the fourteen modular survivors."""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import json

# Phi_60(x)=x^16+x^14-x^10-x^8-x^6+x^2+1.
PHI=[Q(0)]*17
for i,c in {0:1,2:1,6:-1,8:-1,10:-1,14:1,16:1}.items():PHI[i]=Q(c)
def trim(p):
    p=list(p)
    while p and not p[-1]:p.pop()
    return p
def padd(a,b):return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def pmul(a,b):
    p=[Q(0)]*max(0,len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y:p[i+j]+=x*y
    return trim(p)
def divmodp(a,b):
    a=trim(a);b=trim(b);assert b
    q=[Q(0)]*max(0,len(a)-len(b)+1)
    while len(a)>=len(b):
        d=len(a)-len(b);v=a[-1]/b[-1];q[d]=v
        for j,w in enumerate(b):a[j+d]-=v*w
        a=trim(a)
    return trim(q),a
class K:
    def __init__(self,p=0):
        if isinstance(p,K):self.p=p.p;return
        if isinstance(p,(int,Q)):p=[Q(p)]
        self.p=tuple(divmodp(p,PHI)[1])
    def __add__(a,b):return K(padd(a.p,K(b).p))
    __radd__=__add__
    def __neg__(a):return K([-x for x in a.p])
    def __sub__(a,b):return a+-K(b)
    def __rsub__(a,b):return K(b)+-a
    def __mul__(a,b):return K(pmul(a.p,K(b).p))
    __rmul__=__mul__
    def __bool__(a):return bool(a.p)
    def __eq__(a,b):return a.p==K(b).p
    def inverse(a):
        assert a
        r0,r1=PHI,list(a.p);s0,s1=[],[Q(1)]
        while r1:
            q,r=divmodp(r0,r1);r0,r1=r1,r
            s0,s1=s1,padd(s0,[-x for x in pmul(q,s1)])
        assert len(r0)==1
        ans=K([x/r0[0] for x in s0]);assert a*ans==1;return ans
    def __truediv__(a,b):return a*K(b).inverse()
    def conj(a):return sum((c*POW[-j%60] for j,c in enumerate(a.p)),K())
    def serial(a):return [str(x) for x in a.p]
    def approx(a):
        import cmath
        z=cmath.exp(2j*cmath.pi/60)
        return sum(float(x)*z**i for i,x in enumerate(a.p))
POW=[K(1)];z=K([0,1])
for _ in range(59):POW.append(POW[-1]*z)
assert POW[-1]*z==1
for d in (2,3,5):assert POW[60//d]!=1
for i in range(60):assert POW[i]*POW[-i%60]==1

def rref(A,ncols=None):
    a=[[K(x) for x in row] for row in A]
    ncols=len(a[0]) if ncols is None else ncols;piv=[];r=0
    for c in range(ncols):
        s=next((i for i in range(r,len(a)) if a[i][c]),None)
        if s is None:continue
        a[r],a[s]=a[s],a[r];v=a[r][c].inverse();a[r]=[x*v for x in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][c]:
                v=a[i][c];a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        piv.append(c);r+=1
        if r==len(a):break
    return a,piv

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
