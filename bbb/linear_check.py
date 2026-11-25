"""Independent necessary lifting checks by Gaussian elimination over F_7."""
from pathlib import Path
import json

def normal_forms(eq,fix):
    A=[]
    for u,v,w in eq:
        a=[0]*28;a[u]=1;a[v]=-1;a[w]=-1;A.append(a)
    for b in [0,9,18]:
        a=[0]*28
        for j in range(b,b+7):a[j]=1
        A.append(a)
    for i,x in fix:
        a=[0]*28;a[i]=1;a[-1]=x;A.append(a)
    A=[[x%7 for x in a] for a in A]
    piv=[];r=0
    for c in range(27):
        k=next((i for i in range(r,len(A)) if A[i][c]),None)
        if k is None:continue
        A[r],A[k]=A[k],A[r]
        z=pow(A[r][c],-1,7);A[r]=[x*z%7 for x in A[r]]
        for i in range(len(A)):
            if i!=r and A[i][c]:
                z=A[i][c];A[i]=[(x-z*y)%7 for x,y in zip(A[i],A[r])]
        piv.append(c);r+=1
    if any(not any(a[:-1]) and a[-1] for a in A):return None
    free=[j for j in range(27) if j not in piv]
    forms=[None]*27
    for i,j in enumerate(free):
        a=[0]*(len(free)+1);a[i+1]=1;forms[j]=tuple(a)
    for i,j in enumerate(piv):forms[j]=tuple([A[i][-1]]+[-A[i][k]%7 for k in free])
    return forms

def collision(forms):
    if forms is None:return 'inconsistent'
    for b in [0,9,18]:
        for i in range(b,b+7):
            for j in range(i+1,b+7):
                if forms[i]==forms[j]:return [i,j]
    return None

def main():
    root=Path(__file__).resolve().parent
    vals=iter(map(int,(root/'lift_input.txt').read_text().split()))
    n=next(vals);result=[];count=0
    for _ in range(n):
        vi,wi,u0,v0,u1=[next(vals) for _ in range(5)]
        eq=[[next(vals) for _ in range(3)] for _ in range(11)]
        witness=collision(normal_forms(eq,[(u0,0),(v0,0),(u1,1)]))
        result.append({'vi':vi,'wi':wi,'collision':witness})
        count+=witness is not None
    (root/'linear_collisions.json').write_text(json.dumps(result,indent=2)+'\n')
    print(count,'of',n,'excluded by a forced root-label collision.')

if __name__=='__main__':main()
