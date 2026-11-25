"""Two mixed-profile pilots at width 11: common-pentagon CCD and CDF.
Exact arithmetic only. The common-pentagon restriction is explicit:
for CCD all three chosen pentagons coincide; for CDF the C and D
pentagons coincide. These are not full classifications of CCD or CDF.
"""
from pathlib import Path
from itertools import combinations,permutations,product
from fractions import Fraction as F
from collections import Counter
from math import gcd,lcm
import json,time

ROOT=Path(__file__).resolve().parent
U_TRIS=((5,6,7),(8,9,10))
PENT=tuple(range(5))
BASE_C=(PENT,)+U_TRIS

def pairings(a):
    if not a:
        yield ();return
    for b in a[1:]:
        for rest in pairings(tuple(x for x in a[1:] if x!=b)):
            yield ((a[0],b),)+rest

def triples6():
    for rest in combinations(range(6,11),2):
        a=(5,)+rest
        yield (a,tuple(x for x in range(5,11) if x not in a))

def egcd(a,b):
    old_r,r,old_s,s,old_t,t=a,b,1,0,0,1
    while r:
        q=old_r//r
        old_r,r,old_s,s,old_t,t=r,old_r-q*r,s,old_s-q*s,t,old_t-q*t
    if old_r<0:old_r,old_s,old_t=-old_r,-old_s,-old_t
    return old_r,old_s,old_t

def unimodular_reduce(row):
    row=list(row);n=len(row);V=[[int(i==j) for j in range(n)] for i in range(n)]
    for j in range(1,n):
        a,b=row[0],row[j]
        if b==0:continue
        g,s,t=egcd(a,b)
        for v in V:
            x,y=v[0],v[j];v[0],v[j]=s*x+t*y,-(b//g)*x+(a//g)*y
        row[0],row[j]=g,0
    if row[0]<0:
        row[0]=-row[0]
        for v in V:v[0]=-v[0]
    return V,row[0]

def solve_torus(A,b,n):
    states=[(tuple(F(0) for _ in range(n)),tuple(tuple(int(i==j) for j in range(n)) for i in range(n)))]
    for a,target in zip(A,b):
        following=[]
        for c,B in states:
            d=len(B[0]);row=[sum(a[i]*B[i][j] for i in range(n)) for j in range(d)]
            rhs=(target-sum(a[i]*c[i] for i in range(n)))%1
            if not any(row):
                if rhs==0:following.append((c,B))
                continue
            V,g=unimodular_reduce(row)
            BV=[[sum(B[i][k]*V[k][j] for k in range(d)) for j in range(d)] for i in range(n)]
            newB=tuple(tuple(v[1:]) for v in BV)
            for k in range(g):
                newc=tuple((c[i]+BV[i][0]*(rhs+k)/g)%1 for i in range(n))
                following.append((newc,newB))
        states=following
    return states

def row_from_cycles(cycles,phases,orientations,modulus=30):
    row=[0]*11;row[:5]=[modulus*j//5 for j in range(5)]
    for block,phase,sig in zip(cycles,phases,orientations):
        for h,j in enumerate(block):row[j]=(phase+modulus*h*sig//len(block))%modulus
    return row

def cycle_check(exponents,part,modulus):
    for block in part:
        assert sorted((exponents[j]-exponents[block[0]])%modulus for j in block)==list(range(0,modulus,modulus//len(block)))

def three_check(U,V,P,Q,R,M):
    cycle_check(U,P,M);cycle_check(V,Q,M);cycle_check([(a-b)%M for a,b in zip(U,V)],R,M)

def label_permutations():
    return [(0,)+p for p in permutations(range(1,5))]

def ccd():
    # All ten Q partitions, all fifteen R pairings: no support reduction.
    qs=list(triples6());rs=list(pairings(tuple(range(5,11))))
    fifth=[p for p in label_permutations() if len({(j-p[j])%5 for j in range(5)})==5]
    out=[];summary=Counter();seedorders=Counter();feasible=set()
    for qi,q in enumerate(qs):
        uq={j:(k,h) for k,t in enumerate(U_TRIS) for h,j in enumerate(t)}
        vq={j:(k,h) for k,t in enumerate(q) for h,j in enumerate(t)}
        W={}
        for j in range(5,11):
            w=[0]*4;w[uq[j][0]]=1;w[2+vq[j][0]]=-1;W[j]=w
        for ri,r in enumerate(rs):
            A=[[x-y for x,y in zip(W[b],W[a])] for a,b in r]
            for sig in product((-1,1),repeat=4):
                off={j:F(uq[j][1]*sig[uq[j][0]]-vq[j][1]*sig[2+vq[j][0]],3) for j in range(5,11)}
                target=[F(1,2)-off[b]+off[a] for a,b in r]
                for c,B in solve_torus(A,target,4):
                    dim=len(B[0]);M=lcm(30,*(x.denominator for x in c));ph=[int(x*M) for x in c]
                    U=row_from_cycles(U_TRIS,ph[:2],sig[:2],M)
                    V=row_from_cycles(q,ph[2:],sig[2:],M)
                    basisU=[[0]*dim for _ in range(11)];basisV=[[0]*dim for _ in range(11)]
                    for j in range(5,11):
                        basisU[j]=list(B[uq[j][0]]);basisV[j]=list(B[2+vq[j][0]])
                    for pi in fifth:
                        vv=V[:];vv[:5]=[M*x//5 for x in pi]
                        three_check(U,vv,BASE_C,(PENT,)+q,(PENT,)+r,M)
                        diffbasis=[[x-y for x,y in zip(a,b)] for a,b in zip(basisU,basisV)]
                        for bb,part in [(basisU,BASE_C),(basisV,(PENT,)+q),(diffbasis,(PENT,)+r)]:
                            for block in part:assert all(bb[j]==bb[block[0]] for j in block)
                    order=1
                    for row in [U,vv]:
                        for x in row:order=lcm(order,M//gcd(M,x))
                    summary[dim]+=1;seedorders[order]+=1;feasible.add((qi,ri))
                    out.append({"qi":qi,"ri":ri,"Q":q,"R":r,"signs":sig,"modulus":M,
                                "U":U,"V_last6":V[5:],"allowed_fifth_permutations":fifth,
                                "basis_U":basisU,"basis_V":basisV,"dimension":dim,"seed_order":order})
    return {"scope":"all chosen CCD supports sharing a common pentagon",
            "Q_partitions":len(qs),"R_partitions":len(rs),"support_pairs":len(qs)*len(rs),
            "feasible_support_pairs":len(feasible),"phase_components_before_fifth_labels":len(out),
            "dimensions":dict(sorted(summary.items())),"seed_orders":dict(sorted(seedorders.items())),
            "fifth_permutations":fifth,"families":out}

def graph_solver(pairs):
    # Columns outside the pentagon are edges between 2 triangle phases
    # and 3 pair phases. Edge j imposes x_u-y_v=rhs[j].
    pu={j:k for k,t in enumerate(U_TRIS) for j in t}
    pv={j:k+2 for k,t in enumerate(pairs) for j in t}
    edges=[(pu[j],pv[j]) for j in range(5,11)]
    L=[None]*5;L[0]=[0]*6
    while any(x is None for x in L):
        changed=False
        for e,(u,v) in enumerate(edges):
            if L[u] is not None and L[v] is None:
                L[v]=L[u][:];L[v][e]-=1;changed=True
            elif L[v] is not None and L[u] is None:
                L[u]=L[v][:];L[u][e]+=1;changed=True
        assert changed,"disconnected phase graph"
    K=[]
    for e,(u,v) in enumerate(edges):
        row=[x-y for x,y in zip(L[u],L[v])];row[e]-=1
        if any(row):K.append(row)
    assert len(K)==2
    return L,K

def cdf():
    # Q pairings have two orbits under the stabilizer of the two U triangles:
    # one cross-pair, or three cross-pairs.
    qs=(((5,6),(7,8),(9,10)),((5,8),(6,9),(7,10)))
    uinfo={j:(k,h) for k,t in enumerate(U_TRIS) for h,j in enumerate(t)}
    pp=label_permutations()
    rlist=[]
    for anchor in range(5):
        for extra in combinations(range(5,11),2):
            left=[j for j in range(5,11) if j not in extra]
            others=[j for j in range(5) if j!=anchor]
            for perm in permutations(left):
                pair=tuple(zip(others,perm));R=((anchor,)+extra,)+pair
                anchors=[None]*6;er=[0]*6;eps=[0]*6
                for h,j in enumerate(extra,1):anchors[j-5]=anchor;er[j-5]=h
                for a,j in pair:anchors[j-5]=a;eps[j-5]=1
                rlist.append((R,anchors,er,eps))
    assert len(rlist)==1800
    unique={};feasible=set();combos=Counter();tests=0;witnesses=[]
    for qi,q in enumerate(qs):
        L,K=graph_solver(q)
        epsq={j:h for pair in q for h,j in enumerate(pair)}
        for ri,(R,anchors,er,eps) in enumerate(rlist):
            for pi in pp:
                d=[(j-pi[j])%5 for j in range(5)]
                fifth=[6*d[a] for a in anchors]
                for sig in product((-1,1),repeat=3):
                    rhs=[(fifth[j-5]+10*(er[j-5]*sig[2]-uinfo[j][1]*sig[uinfo[j][0]])+
                          15*(eps[j-5]+epsq[j]))%30 for j in range(5,11)]
                    tests+=1
                    if any(sum(a*b for a,b in zip(k,rhs))%30 for k in K):continue
                    theta=[sum(a*b for a,b in zip(l,rhs))%30 for l in L]
                    U=row_from_cycles(U_TRIS,theta[:2],sig[:2])
                    V=row_from_cycles(q,theta[2:],(1,1,1));V[:5]=[6*x for x in pi]
                    three_check(U,V,BASE_C,(PENT,)+q,R,30)
                    # For a generic common phase on the last six columns,
                    # alternative C/D decompositions must respect the blocks.
                    u_also_D=all((x+15)%30 in U[5:] for x in U[5:]) and len(set(U[5:]))==6
                    v_also_C=len(set(V[5:]))==6 and set((x-V[5])%30 for x in V[5:])==set(range(0,30,5))
                    category=f"U_also_D={u_also_D},V_also_C={v_also_C}"
                    combos[category]+=1;feasible.add((qi,ri))
                    key=tuple(U+V);record={"qi":qi,"ri":ri,"Q":q,"R":R,"pi":pi,"signs":sig,
                            "modulus":30,"U":U,"V":V,"dimension":1,
                            "U_also_D":u_also_D,"V_also_C":v_also_C}
                    unique.setdefault(key,record)
        print("CDF Q orbit",qi+1,"feasible witnesses",sum(combos.values()),flush=True)
    return {"scope":"all chosen CDF supports with coincident C and D pentagons",
            "Q_pairing_orbits":2,"R_partitions_per_orbit":1800,"support_presentations":3600,
            "fifth_label_permutations":len(pp),"orientation_choices":8,
            "tested_combinations":tests,"feasible_support_presentations":len(feasible),
            "feasible_label_orientation_choices":sum(combos.values()),
            "overlap_categories":dict(combos),"distinct_normalized_seeds":len(unique),
            "families":list(unique.values())}

def main():
    t=time.time();a=ccd()
    (ROOT/"ccd_common_pentagon.json").write_text(json.dumps(a,indent=2)+"\n")
    print("CCD",{k:v for k,v in a.items() if k!="families"},flush=True)
    b=cdf()
    (ROOT/"cdf_common_pentagon.json").write_text(json.dumps(b,indent=2)+"\n")
    print("CDF",{k:v for k,v in b.items() if k!="families"},flush=True)
    print("seconds",round(time.time()-t,3),flush=True)

if __name__=="__main__":
    main()

