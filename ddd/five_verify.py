"""Independent quintic-label verification: Gaussian elimination and enumeration.

Does not call or reproduce the C++ domain-propagation search.
"""
from itertools import product

def forms(eq):
    A=[]
    for u,v,w in eq:
        a=[0]*25;a[u]=1;a[v]=-1;a[w]=-1;A.append(a)
    for b in (0,8,16):
        a=[0]*25
        for j in range(b,b+5):a[j]=1
        A.append(a)
    for i,x in ((0,0),(8,0),(1,1)):
        a=[0]*25;a[i]=1;a[-1]=x;A.append(a)
    A=[[x%5 for x in a] for a in A];piv=[];r=0
    for c in range(24):
        k=next((i for i in range(r,len(A)) if A[i][c]),None)
        if k is None:continue
        A[r],A[k]=A[k],A[r];z=pow(A[r][c],-1,5)
        A[r]=[x*z%5 for x in A[r]]
        for i in range(len(A)):
            if i!=r and A[i][c]:
                z=A[i][c];A[i]=[(x-z*y)%5 for x,y in zip(A[i],A[r])]
        piv.append(c);r+=1
    if any(not any(a[:-1]) and a[-1] for a in A):return None
    free=[j for j in range(24) if j not in piv];out=[None]*24
    for i,j in enumerate(free):
        a=[0]*(len(free)+1);a[i+1]=1;out[j]=tuple(a)
    for i,j in enumerate(piv):out[j]=tuple([A[i][-1]]+[-A[i][k]%5 for k in free])
    return out

def collision(F):
    if F is None:return [-1,-1]
    for b in (0,8,16):
        for i in range(b,b+5):
            for j in range(i+1,b+5):
                if F[i]==F[j]:return [i,j]
    return None

def exhaustive(F):
    """Return a witness or None by exhausting the affine linear solution space."""
    if F is None:return None,0
    d=len(F[0])-1;tested=0
    for x in product(range(5),repeat=d):
        tested+=1
        values=[(a[0]+sum(c*t for c,t in zip(a[1:],x)))%5 for a in F]
        if all(len(set(values[b:b+5]))==5 for b in (0,8,16)):
            return values,tested
    return None,tested
