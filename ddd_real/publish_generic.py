"""Recover the exact printed order: tail factor 1, a, b, then coefficients."""
from pathlib import Path
import json
from verify_rows import verify
R=Path(__file__).resolve().parent
def main():
 rows=[tuple(map(int,l.split()))for l in(R/'analysis/generic_records.txt').read_text().splitlines()]
 factors={(0,0):0,(1,0):1,(0,1):2}
 rows.sort(key=lambda z:(factors[z[48:50]],tuple(z[16::5])))
 assert len(rows)==24
 path=R/'analysis/generic_published.txt'
 path.write_text(''.join(' '.join(map(str,z))+'\n'for z in rows))
 result=verify(path,R/'b_equals_c.txt')
 (R/'analysis/generic_published.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
