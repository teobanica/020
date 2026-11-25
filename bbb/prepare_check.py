from pathlib import Path
import json
from collapse import partitions,BASE,equations
root=Path(__file__).resolve().parent
d=json.loads((root/'collapsed.json').read_text());parts=partitions()
lines=[];count=0
for vi,o in enumerate(d['support_orbits']):
    for wi,R in enumerate(parts):
        if not set(BASE[0])&set(o['representative'][0])&set(R[0]):continue
        A,b=equations(o['representative'],R)
        lines.append(f'{vi} {wi} {len(A)}')
        for ar,rhs in zip(A,b):
            assert (2*rhs).denominator==1
            lines.append(' '.join(map(str,[*ar,int(2*rhs)])))
        count+=1
(root/'check_input.txt').write_text('\n'.join([str(count)]+lines)+'\n')
