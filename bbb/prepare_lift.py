from pathlib import Path
import json
from collapse import partitions, BASE
from linear_check import normal_forms,collision

ROOT=Path(__file__).resolve().parent

def column_variables(P,base):
    a=[None]*11
    for k,col in enumerate(P[0]):a[col]=base+k
    for k,pair in enumerate(P[1:]):
        for col in pair:a[col]=base+7+k
    return a

def main():
    data=json.loads((ROOT/'collapsed.json').read_text());parts=partitions()
    pairs=sorted({(r['vi'],r['wi']) for r in data['records']})
    out=[];kept=0;excluded=[]
    for vi,wi in pairs:
        P=data['support_orbits'][vi]['representative'];R=parts[wi]
        u=column_variables(BASE,0);v=column_variables(P,9);w=column_variables(R,18)
        col=min(set(BASE[0])&set(P[0])&set(R[0]))
        other=next(j for j in BASE[0] if j!=col)
        witness=collision(normal_forms(list(zip(u,v,w)),[(u[col],0),(v[col],0),(u[other],1)]))
        if witness is not None:
            excluded.append({'vi':vi,'wi':wi,'collision':witness});continue
        kept+=1
        out.append(' '.join(map(str,(vi,wi,u[col],v[col],u[other]))))
        out.extend(' '.join(map(str,t)) for t in zip(u,v,w))
    (ROOT/'lift_input.txt').write_text('\n'.join([str(kept)]+out)+'\n')
    (ROOT/'excluded_lifts.json').write_text(json.dumps(excluded,indent=2)+'\n')
    print('Prepared',kept,'prime-7 lifting problems;',len(excluded),'excluded linearly.')

if __name__=='__main__':main()
