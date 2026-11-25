"""Exact independent certificate verification and symbolic cycle verification."""
from pathlib import Path
from itertools import product
from math import prod,gcd,lcm
from collections import Counter
import json,time
from supports import BASE,partitions,equations,row,orbits
from prepare_five import variables
from prepare_three import orientation_matrix,offsets
from five_verify import forms,collision,exhaustive
ROOT=Path(__file__).resolve().parent

def mm(A,B):
    cols=list(zip(*B))
    return [[sum(x*y for x,y in zip(a,b)) for b in cols] for a in A]

def det(A):
    a=[r[:] for r in A];n=len(a);old=1;sign=1
    for k in range(n-1):
        if not a[k][k]:
            j=next((i for i in range(k+1,n) if a[i][k]),None)
            if j is None:return 0
            a[k],a[j]=a[j],a[k];sign=-sign
        z=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=z*a[i][j]-a[i][k]*a[k][j]
                assert num%old==0
                a[i][j]=num//old
        for i in range(k+1,n):a[i][k]=0
        old=z
    return sign*a[-1][-1]

def exact_family(Q,R,s,mask,branch,labels):
    V=s['V'];rank=s['rank'];diag=s['diagonal'];U=s['U']
    sig=[1 if mask>>j&1 else -1 for j in range(6)]
    coeff=orientation_matrix(Q,R)
    rhs=[sum(a*b for a,b in zip(r,sig)) for r in coeff]
    ub=[sum(a*b for a,b in zip(r,rhs)) for r in U]
    theta=[sum(V[i][j]*(30//diag[j])*(ub[j]+3*branch[j]) for j in range(rank))%90 for i in range(4)]
    seed=[];basis=[]
    for P,off,base in ((BASE,0,0),(Q,2,7)):
        D=row(P,off);eps=offsets(P,off);lab=variables(P,5,base)
        seed.append([(18*labels[lab[j]]+30*sum(a*b for a,b in zip(eps[j],sig))+sum(a*b for a,b in zip(D[j],theta)))%90 for j in range(11)])
        basis.append([[sum(D[j][i]*V[i][k] for i in range(4)) for k in range(rank,4)] for j in range(11)])
    ws=[(a-b)%90 for a,b in zip(*seed)]
    wb=[[a-b for a,b in zip(x,y)] for x,y in zip(*basis)]
    for const,coeff,P in zip(seed+[ws],basis+[wb],(BASE,Q,R)):
        for block in P:
            j=block[0];p=len(block)
            assert sorted((const[k]-const[j])%90 for k in block)==[90*k//p for k in range(p)]
            assert all(coeff[k]==coeff[j] for k in block)
    return {'modulus':90,'exponents':seed,'basis':basis,'supports':[BASE,Q,R]}

def main():
    start=time.time();data=json.loads((ROOT/'supports.json').read_text());parts=partitions();reps=orbits(parts)
    assert len(parts)==4620 and len(set(parts))==4620 and len(reps)==21
    assert [(tuple(tuple(b) for b in o['representative']),o['size']) for o in data['support_orbits']]==reps
    five=json.loads((ROOT/'five.json').read_text());phase=json.loads((ROOT/'phase_check.json').read_text());three=json.loads((ROOT/'three.json').read_text())
    assert len(five)==len(phase)==len(three)==97020
    stats=Counter();dimensions=Counter();supportdims=Counter();orders=Counter();examples={};full_search_points=0
    for idx,(f,p,s) in enumerate(zip(five,phase,three)):
        vi,wi=divmod(idx,4620);assert all((r['vi'],r['wi'])==(vi,wi) for r in (f,p,s))
        Q=data['support_orbits'][vi]['representative'];R=parts[wi]
        A=equations(Q,R);U=s['U'];V=s['V'];rank=s['rank'];diag=s['diagonal']
        assert abs(det(U))==abs(det(V))==1
        D=mm(mm(U,A),V)
        assert D==[[diag[i] if i==j and i<rank else 0 for j in range(4)] for i in range(8)]
        assert all(x>0 for x in diag) and all(b%a==0 for a,b in zip(diag,diag[1:]))
        assert rank==p['rank'] and prod(diag)==p['torsion_order'] and (diag[-1] if diag else 1)==p['max_smith_invariant']
        assert all(x in (1,2,3) for x in diag) and prod(diag) in (1,2,3)
        # All possible triangle orientations, independently re-evaluated.
        UB=mm(U,orientation_matrix(Q,R));allowed=[]
        for mask in range(64):
            sig=[1 if mask>>j&1 else -1 for j in range(6)]
            if all(sum(a*b for a,b in zip(rr,sig))%3==0 for rr in UB[rank:]):allowed.append(mask)
        assert allowed==s['orientations']
        maps=[variables(P,5,b) for P,b in ((BASE,0),(Q,7),(R,14))];eq=list(zip(*maps))
        if f['satisfiable']:
            labels=f['labels'];assert len(labels)==21 and all(0<=x<5 for x in labels)
            assert all(sorted(labels[b:b+5])==list(range(5)) for b in (0,7,14))
            assert all((labels[a]-labels[b]-labels[c])%5==0 for a,b,c in eq)
            stats['five_feasible']+=1
        else:
            F=forms(eq)
            if 'collision' in f:
                assert collision(F)==f['collision'];stats['five_linear_excluded']+=1
            else:
                witness,n=exhaustive(F);assert witness is None;full_search_points+=n
                stats['five_exhaustively_excluded']+=1
            continue
        if not allowed:stats['three_excluded_after_five']+=1;continue
        stats['feasible_supports']+=1;supportdims[4-rank]+=1;orders[prod(diag)]+=1
        for mask in allowed:
            for branch in product(*(range(d) for d in diag)):
                F=exact_family(Q,R,s,mask,branch,labels)
                dimensions[4-rank]+=1
                # Cross-ratio exponents include row ratios relative to column 0.
                cross_order=lcm(*(90//gcd(90,(r[j]-r[0])%90) for r in F['exponents'] for j in range(11)))
                if cross_order not in examples:
                    g=90//cross_order
                    ex=[[(x-r[0])%90//g for x in r] for r in F['exponents']]
                    examples[cross_order]={'vi':vi,'wi':wi,'mask':mask,'branch':branch,'dimension':4-rank,
                        'modulus':cross_order,'exponents':ex,'supports':F['supports'],'basis':F['basis']}
        if idx%5000==0:print('Verified',idx,'systems',flush=True)
    result={'verified':True,'partitions':len(parts),'relative_support_orbits':len(reps),'support_pairs':len(five),
        **dict(stats),'exhaustive_quintic_points':full_search_points,'feasible_supports_by_dimension':dict(sorted(supportdims.items())),
        'phase_components_by_dimension':dict(sorted(dimensions.items())),
        'phase_components':sum(dimensions.values()),'feasible_supports_by_torsion_order':dict(sorted(orders.items())),
        'example_root_orders':sorted(examples),'note':'Counts are redundant support/orientation/phase presentations, not equivalence classes. One quintic witness per feasible support; all quintic labels are permitted by the normal form.',
        'seconds':round(time.time()-start,3)}
    (ROOT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'examples.json').write_text(json.dumps(examples,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
