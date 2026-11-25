from pathlib import Path
import json
from supports import BASE,partitions,equations
ROOT=Path(__file__).resolve().parent
def offsets(P,start):
    a=[[0]*6 for _ in range(11)]
    for k,block in enumerate(P[1:]):
        for j,col in enumerate(block):a[col][start+k]=j
    return a
def orientation_matrix(Q,R):
    p,q,r=offsets(BASE,0),offsets(Q,2),offsets(R,4)
    z=[[r[j][k]-p[j][k]+q[j][k] for k in range(6)] for j in range(11)]
    return [[a-b for a,b in zip(z[j],z[block[0]])] for block in R for j in block[1:]]
def main():
    data=json.loads((ROOT/'supports.json').read_text());parts=partitions()
    keep=[{'vi':vi,'wi':wi} for vi in range(len(data['support_orbits'])) for wi in range(len(parts))]
    with (ROOT/'three_input.txt').open('w') as f:
        f.write(str(len(keep))+'\n')
        for r in keep:
            Q=data['support_orbits'][r['vi']]['representative'];R=parts[r['wi']]
            f.write(f"{r['vi']} {r['wi']}\n")
            for a,b in zip(equations(Q,R),orientation_matrix(Q,R)):
                f.write(' '.join(map(str,a+b))+'\n')
    print('Prepared',len(keep),'ternary phase systems.',flush=True)
if __name__=='__main__':main()
