from pathlib import Path
import sys,json,itertools
root=Path(__file__).resolve().parent;old=root.parent/'n11_noA_DDD';sys.path.insert(0,str(old))
from supports import BASE,partitions,equations,row,offsets
classes=json.loads((old/'support_classes.json').read_text())['classes'];phases=json.loads((old/'phases.json').read_text());parts=partitions();orbits=json.loads((old/'supports.json').read_text())['support_orbits']
records=[]
for idx,c in enumerate(classes):
 if not c['feasible']:continue
 vi,wi=c['representative'];Q=orbits[vi]['representative'];R=parts[wi];ph=phases[vi*6930+wi];assert(ph['vi'],ph['wi'])==(vi,wi);assert ph['binary_feasible'];A,g=equations(Q,R)
 records.append({'class_id':idx,'vi':vi,'wi':wi,'supports':[BASE,Q,R],'A':A,'g':g,'rank':ph['rank'],'dimension':c['dimension'],'diagonal':ph['diagonal'],'U':ph['U'],'V':ph['V']})
(root/'classes.json').write_text(json.dumps(records,indent=2)+'\n')
with(root/'labels_input.txt').open('w')as f:
 f.write(str(len(records))+'\n')
 for r in records:
  f.write(' '.join(map(str,[r['class_id'],r['dimension'],r['rank']]+r['diagonal']))+'\n')
  for s in r['supports']:f.write(' '.join(map(str,[j for b in s for j in b]))+'\n')
  for name in ('U','V'):
   for line in r[name]:f.write(' '.join(map(str,line))+'\n')
  f.write(' '.join(map(str,r['g']))+'\n')
print('classes',len(records));print('exception',[r['class_id']for r in records if 4 in r['diagonal']]);print('dim3',[(r['class_id'],r['supports'])for r in records if r['dimension']==3])
