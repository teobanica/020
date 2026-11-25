"""Independent checks for the two common-pentagon pilots.
CCD: exhaustive sixth-root phase grid and primitive integer family bases.
CDF: direct antipodal identities, independent of the graph solver.
Profile audit: exhaustive exact covers by prime cycles.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,permutations,product
from collections import Counter
from math import gcd
import json
ROOT=Path(__file__).resolve().parent
PROFILES={"A":(11,),"B":(7,2,2),"C":(5,3,3),"D":(5,2,2,2),
          "E":(3,3,3,2),"F":(3,2,2,2,2)}

def determinant(a):
    n=len(a)
    if n==0:return 1
    if n==1:return a[0][0]
    return sum((-1)**j*a[0][j]*determinant([row[:j]+row[j+1:] for row in a[1:]]) for j in range(n))

def rank(a):
    a=[[F(x) for x in row] for row in a];r=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][j]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]
        z=a[r][j];a[r]=[x/z for x in a[r]]
        for i in range(r+1,len(a)):
            z=a[i][j];a[i]=[x-z*y for x,y in zip(a[i],a[r])]
        r+=1
    return r

def tri_partitions():
    for t in combinations(range(6,11),2):
        a=(5,)+t;yield (a,tuple(x for x in range(5,11) if x not in a))

def pair_partitions(a):
    if not a:yield ();return
    for j in range(1,len(a)):
        for r in pair_partitions(a[1:j]+a[j+1:]):
            yield ((a[0],a[j]),)+r

def verify_ccd():
    data=json.loads((ROOT/"ccd_common_pentagon.json").read_text())
    qs=list(tri_partitions());rs=list(pair_partitions(tuple(range(5,11))))
    families={}
    for rec in data["families"]:
        k=(rec["qi"],rec["ri"],tuple(rec["signs"]))
        families.setdefault(k,[]).append(rec)
    total=0;dimensions=Counter();supports=set();gridpoints=0;visited=0
    grids=list(product(range(6),repeat=4))
    Utri=((5,6,7),(8,9,10))
    for qi,q in enumerate(qs):
        uinfo={j:(b,h) for b,t in enumerate(Utri) for h,j in enumerate(t)}
        vinfo={j:(b,h) for b,t in enumerate(q) for h,j in enumerate(t)}
        vectors={}
        for j in range(5,11):
            x=[0]*4;x[uinfo[j][0]]=1;x[2+vinfo[j][0]]=-1;vectors[j]=x
        for ri,r in enumerate(rs):
            A=[[x-y for x,y in zip(vectors[b],vectors[a])] for a,b in r]
            rr=rank(A)
            for sig in product((-1,1),repeat=4):
                key=(qi,ri,sig);actual=set()
                constants={j:2*(sig[uinfo[j][0]]*uinfo[j][1]-sig[2+vinfo[j][0]]*vinfo[j][1]) for j in range(5,11)}
                rhs=[(3-constants[b]+constants[a])%6 for a,b in r]
                for x in grids:
                    visited+=1
                    if all(sum(u*v for u,v in zip(a,x))%6==target for a,target in zip(A,rhs)):
                        actual.add(x)
                covered=set()
                for rec in families.get(key,[]):
                    dim=rec["dimension"];assert dim==4-rr
                    seed=[rec["U"][t[0]]*6//rec["modulus"] for t in Utri]
                    seed+=[rec["V_last6"][t[0]-5]*6//rec["modulus"] for t in q]
                    assert all(rec["U"][t[0]]*6%rec["modulus"]==0 for t in Utri)
                    assert all(rec["V_last6"][t[0]-5]*6%rec["modulus"]==0 for t in q)
                    B=[rec["basis_U"][t[0]] for t in Utri]+[rec["basis_V"][t[0]] for t in q]
                    minor_gcd=0
                    for rows in combinations(range(4),dim):
                        minor_gcd=gcd(minor_gcd,abs(determinant([B[i] for i in rows])))
                    assert minor_gcd==1
                    assert all(sum(a[i]*B[i][j] for i in range(4))==0 for a in A for j in range(dim))
                    pts={tuple((seed[i]+sum(B[i][j]*t[j] for j in range(dim)))%6 for i in range(4))
                         for t in product(range(6),repeat=dim)}
                    assert len(pts)==6**dim and not (covered&pts)
                    covered|=pts;dimensions[dim]+=1;total+=1;supports.add((qi,ri))
                assert covered==actual
                gridpoints+=len(actual)
    # Rank and all primitive components are independently checked. Completeness
    # on the continuous torus also uses the exhaustive unimodular branch proof
    # in pilots.py; a finite grid alone is not asserted to prove completeness.
    assert visited==150*16*6**4 and total==684
    assert dimensions==Counter({2:648,3:36}) and len(supports)==96
    fifth=[]
    for p in permutations(range(1,5)):
        pi=(0,)+p
        if sorted((j-pi[j])%5 for j in range(5))==list(range(5)):fifth.append(pi)
    assert [list(x) for x in fifth]==data["fifth_permutations"]
    return {"verified":True,"grid_points_tested":visited,"grid_solutions":gridpoints,
            "phase_components":total,"feasible_supports":len(supports),
            "dimensions":dict(sorted(dimensions.items())),"primitive_bases":True,
            "fifth_label_assignments":len(fifth)}

def all_profiles(row,M=30,basis=None):
    if basis is None:basis=[()] * len(row)
    cycles={}
    for p in (2,3,5,7,11):
        cycles[p]=[]
        if M%p:continue
        for c in combinations(range(len(row)),p):
            if not all(basis[j]==basis[c[0]] for j in c):continue
            if sorted((row[j]-row[c[0]])%M for j in c)==list(range(0,M,M//p)):
                cycles[p].append((sum(1<<j for j in c),list(c)))
    out={}
    for name,profile in PROFILES.items():
        def rec(k,used,blocks):
            if k==len(profile):
                return blocks if used==(1<<len(row))-1 else None
            for mask,c in cycles[profile[k]]:
                if used&mask:continue
                z=rec(k+1,used|mask,blocks+[c])
                if z is not None:return z
            return None
        z=rec(0,0,[])
        if z is not None:out[name]=z
    return out

def verify_cdf():
    data=json.loads((ROOT/"cdf_common_pentagon.json").read_text())
    expected={tuple(r["U"]+r["V"]) for r in data["families"]}
    found=set();count=0;feasible=set();byorbit=Counter()
    utri=((5,6,7),(8,9,10));ui={j:(k,h) for k,t in enumerate(utri) for h,j in enumerate(t)}
    qs=(((5,6),(7,8),(9,10)),((5,8),(6,9),(7,10)))
    pis=[(0,)+p for p in permutations(range(1,5))]
    tested=0;strict=None
    for qi,q in enumerate(qs):
        ri=-1
        for anchor in range(5):
            for extra in combinations(range(5,11),2):
                left=[j for j in range(5,11) if j not in extra]
                others=[j for j in range(5) if j!=anchor]
                for perm in permutations(left):
                    ri+=1
                    for pi in pis:
                        fifth=[6*(j-pi[j])%30 for j in range(5)]
                        for su,sv,sw in product((-1,1),repeat=3):
                            tested+=1
                            w=fifth+[0]*6
                            w[extra[0]]=(fifth[anchor]+10*sw)%30
                            w[extra[1]]=(fifth[anchor]+20*sw)%30
                            for a,j in zip(others,perm):w[j]=(fifth[a]+15)%30
                            signs=(su,sv);e={j:10*signs[ui[j][0]]*ui[j][1] for j in range(5,11)}
                            delta=None;ok=True
                            for a,b in q:
                                if ui[a][0]==ui[b][0]:
                                    if (e[b]-e[a]-w[b]+w[a])%30!=15:ok=False;break
                                else:
                                    assert ui[a][0]==0 and ui[b][0]==1
                                    target=(15+w[b]-w[a]+e[a]-e[b])%30
                                    if delta is not None and target!=delta:ok=False;break
                                    delta=target
                            if not ok:continue
                            assert delta is not None
                            U=[6*j for j in range(5)]+[(e[j]+(delta if ui[j][0] else 0))%30 for j in range(5,11)]
                            V=[6*j for j in pi]+[(U[j]-w[j])%30 for j in range(5,11)]
                            found.add(tuple(U+V));count+=1;feasible.add((qi,ri));byorbit[qi]+=1
    assert tested==691200 and count==data["feasible_label_orientation_choices"]==1920
    assert len(feasible)==data["feasible_support_presentations"]==480
    assert found==expected and len(found)==168
    candidate=next(r for r in data["families"] if not r["U_also_D"] and not r["V_also_C"])
    U,V=candidate["U"],candidate["V"];w=[(x-y)%30 for x,y in zip(U,V)]
    actual=[all_profiles(x) for x in (U,V,w)]
    assert [set(x) for x in actual]==[{"C"},{"D"},{"F"}]
    generic_basis=[(0,)]*5+[(1,)]*6
    generic=[all_profiles(U,basis=generic_basis),all_profiles(V,basis=generic_basis),all_profiles(w)]
    assert [set(x) for x in generic]==[{"C"},{"D"},{"F"}]
    # CCD explicit family: rows eta^j and eta^(2j), then s*T,s*a*T
    # and t*T,-t*a*T. Take s=t=1 and a=eta for a strict root-30 point.
    CU=[0,6,12,18,24,0,10,20,6,16,26]
    CV=[0,12,24,6,18,0,10,20,21,1,11]
    CW=[(x-y)%30 for x,y in zip(CU,CV)]
    ca=[all_profiles(x) for x in (CU,CV,CW)]
    assert [set(x) for x in ca]==[{"C"},{"C"},{"D"}]
    z=(0,0,0);bu=[z]*5+[(1,0,0)]*3+[(1,0,1)]*3
    bv=[z]*5+[(0,1,0)]*3+[(0,1,1)]*3
    # At the phase-zero seed a=s=t=1 the second triangle constants change;
    # here a=eta is absorbed in the constant, leaving the same parameter basis.
    bw=[tuple(x-y for x,y in zip(a,b)) for a,b in zip(bu,bv)]
    cg=[all_profiles(CU,basis=bu),all_profiles(CV,basis=bv),all_profiles(CW,basis=bw)]
    assert [set(x) for x in cg]==[{"C"},{"C"},{"D"}]
    report={"verified":True,"tested_combinations":tested,"feasible_choices":count,
            "by_Q_orbit":dict(byorbit),"feasible_supports":len(feasible),"normalized_seeds":len(found),
            "strict_CDF_example":candidate,"strict_CDF_seed_profiles":actual,
            "strict_CDF_generic_profiles":generic,"strict_CCD_example":{"modulus":30,"U":CU,"V":CV,
            "basis_U":bu,"basis_V":bv,"seed_profiles":ca,"generic_profiles":cg}}
    return report

def verify_five_rows():
    # Parameters (a,b,s,t,r,z). Constants are exponents of a 30th root.
    zero=(0,)*6
    rows=[];basis=[]
    for i,(mode,sign,amp) in enumerate([(0,1,None),(1,1,0),(1,-1,0),(2,1,1),(2,-1,1)]):
        e=[6*i*j%30 for j in range(5)];bb=[zero]*5
        phase=None if i==0 else i+1
        for block in range(2):
            for j in range(3):
                powers=[0]*6
                if phase is not None:powers[phase]=1
                if block and amp is not None:powers[amp]=1
                e.append((10*mode*j+(15 if block and sign==-1 else 0))%30)
                bb.append(tuple(powers))
        rows.append(e);basis.append(bb)
    types={};certificates=[]
    for i,j in combinations(range(5),2):
        d=[(a-b)%30 for a,b in zip(rows[i],rows[j])]
        b=[tuple(x-y for x,y in zip(a,c)) for a,c in zip(basis[i],basis[j])]
        profiles=all_profiles(d,basis=b)
        expected="D" if (i,j) in [(1,2),(3,4)] else "C"
        assert set(profiles)=={expected}
        types[i,j]=expected;certificates.append({"rows":[i,j],"profile":expected,"cycles":profiles[expected]})
    tri=Counter()
    for a,b,c in combinations(range(5),3):
        tri["".join(sorted([types[a,b],types[a,c],types[b,c]]))]+=1
    assert tri==Counter({"CCD":6,"CCC":4})
    # Recover all six parameters from dephased entries and ratios.
    vectors=[basis[i][5] for i in range(1,5)]
    vectors += [tuple(x-y for x,y in zip(basis[i][8],basis[i][5])) for i in [1,3]]
    assert set(vectors)=={tuple(int(i==j) for i in range(6)) for j in range(6)}
    result={"verified":True,"dimension":6,"modulus":30,"exponents":rows,"basis":basis,
            "quotient_certificates":certificates,"generic_three_row_types":dict(tri),
            "note":"Every member extends the displayed three-parameter CCD family; generic edge types C except (1,2) and (3,4), which are D."}
    (ROOT/"five_row_family.json").write_text(json.dumps(result,indent=2)+"\n")
    return {"verified":True,"dimension":6,"generic_triple_types":dict(tri)}

def main():
    a=verify_ccd();print("CCD independent verification",a,flush=True)
    b=verify_cdf();print("CDF independent verification",{k:v for k,v in b.items() if "example" not in k and "profiles" not in k},flush=True)
    c=verify_five_rows();print("Five-row extension",c,flush=True)
    out={"verified":True,"CCD":a,"CDF":b,"five_row_CCD_extension":c,
         "scope":"complete common-pentagon pilots only; explicit CCD extension to five rows; other pentagon supports and CDF extensions are untested"}
    (ROOT/"verification.json").write_text(json.dumps(out,indent=2)+"\n")

if __name__=="__main__":main()

