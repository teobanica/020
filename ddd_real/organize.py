from pathlib import Path
from collections import defaultdict,Counter
from math import gcd,lcm
import json
base=Path(__file__).resolve().parent;source=base/'exceptional_equations_previous.txt'
components=defaultdict(set);equations=[]
for line in source.read_text().splitlines():
 r,s,t,e=map(int,line.split());d=gcd(gcd(abs(r),abs(s)),abs(t));v=(r//d,s//d,t//d);assert e%2==0
 # Every component lies in mu_840 in this particular saved list.
 roots=[f for f in range(840)if (d*f-e)%840==0]
 assert len(roots)==d,(r,s,t,e,d)
 for f in roots:components[v].add(f)
 equations.append((r,s,t,e))
report={'equations':len(equations),'directions':len(components),'connected_tori':sum(map(len,components.values())),'tori':[{'normal':list(v),'roots_mod840':sorted(es)}for v,es in sorted(components.items())]}
(base/'connected_tori.json').write_text(json.dumps(report,indent=2)+'\n')
print('equations',len(equations),'directions',len(components),'connected tori',report['connected_tori'])
print([(v,len(es),dict(Counter(840//gcd(840,e)for e in es)))for v,es in sorted(components.items())])
