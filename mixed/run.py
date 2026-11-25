"""Reproduce both restricted pilots and all exact checks."""
from pathlib import Path
import subprocess,sys,json

ROOT=Path(__file__).resolve().parent

if __name__=="__main__":
    for filename in ["pilots.py","verify.py"]:
        subprocess.run([sys.executable,str(ROOT/filename)],cwd=ROOT,check=True)
    result=json.loads((ROOT/"verification.json").read_text())
    assert result["verified"]
    assert result["CCD"]["phase_components"]==684
    assert result["CDF"]["normalized_seeds"]==168
    assert result["five_row_CCD_extension"]["dimension"]==6
    print("All CCD/CDF pilot checks passed.")
