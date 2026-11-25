"""Reproduce the FFF census and its checks. Standard library only."""
from pathlib import Path
from collections import Counter
import argparse,json,subprocess,sys
ROOT=Path(__file__).resolve().parent

def read(name):
    return json.loads((ROOT/name).read_text())

def call(args):
    print("+"," ".join(map(str,args)),flush=True)
    subprocess.run(list(map(str,args)),cwd=ROOT,check=True)

def prepare_ring_input():
    s=read("fff_supports.json");f=read("fff_feasible.json")
    with (ROOT/"verification_input.txt").open("w") as out:
        print(len(s["orbits"]),len(s["partitions"]),len(f),file=out)
        for r in s["orbits"]:
            print(*sum(r["representative"],[]),file=out)
        for r in s["partitions"]:
            print(*sum(r,[]),file=out)
        for row in f:
            print(*row,file=out)

def prepare_class_input():
    c=read("fff_classes.json")["classes"]
    with (ROOT/"class_verify_input.txt").open("w") as out:
        for r in c:
            print(r["tensor"],*(r[k] for k in
                  ["dimension","largest_smith","torsion","orientations",
                   "presentations","raw_tensors"]),file=out)

def summarize():
    s=read("fff_summary.json");c=read("fff_classes.json")
    r=read("ring_verification.json");g=read("class_verification.json")
    f=read("family_verification.json");e=read("example_verification.json")
    assert s["partitions"]==17325 and s["relative_orbits"]==44
    assert s["presentations"]==r["presentations"]==762300
    assert s["binary"]==607052
    assert s["feasible"]==r["feasible_presentations"]==g["presentations"]==552102
    assert s["orientations"]==r["orientations"]==3876166
    assert s["components"]==r["phase_components"]==4933076
    assert s["raw_tensors"]==g["raw_tensors"]==320290
    assert c["feasible_support_classes"]==g["classes"]==f["support_classes"]==19841
    assert r["verified"] and g["verified"] and f["verified"] and e["verified"]
    dims={str(k):v for k,v in sorted(Counter(x["dimension"] for x in c["classes"]).items())}
    assert dims==g["dimensions"]=={"2":17253,"3":2394,"4":175,"5":17,"6":2}
    assert f["orientations"]==146608 and f["components"]==190328
    assert f["seed_orders"]=={"6":173112,"12":15040,"24":2176}
    assert s["orders"]=={"6":4550428,"12":347352,"24":35296}
    assert s["integer_audit"]["matrices"]==762300
    assert sum(x["presentations"] for x in c["classes"])==552102
    assert sum(x["raw_tensors"] for x in c["classes"])==320290
    result={"verified":True,"profile":[3,2,2,2,2],
            "support_presentations":762300,"feasible_presentations":552102,
            "support_classes":19841,"classes_by_dimension":dims,
            "affine_presentations_for_class_representatives":190328,
            "seed_orders":[6,12,24],
            "checks":["Smith identities and unimodularity for every presentation",
                      "independent finite-ring feasibility for every presentation and orientation",
                      "full-group canonical minima and transformation witnesses",
                      "symbolic cycle and affine-basis check of every retained component",
                      "explicit six-parameter family and primitive-order-24 example"],
            "scope":"chosen-support normal form; not full matrix equivalence or row extensions"}
    (ROOT/"verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    return result

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--check-existing",action="store_true",
                        help="cross-check the supplied executed verification outputs")
    args=parser.parse_args()
    if not args.check_existing:
        call(["node","fff.js"])
        for name in ["classes","verify_rings","verify_classes"]:
            call(["g++","-std=c++17","-O3",name+".cpp","-o",name])
        prepare_ring_input()
        call(["./verify_rings","verification_input.txt","ring_verification.json"])
        call(["./classes","class_input.txt","fff_classes.json","class_witnesses.txt"])
        prepare_class_input()
        call(["./verify_classes","class_verify_input.txt","class_input.txt",
              "class_witnesses.txt","class_verification.json"])
        call(["node","verify_families.js"])
        call([sys.executable,"verify_examples.py"])
    summarize()

if __name__=="__main__":
    main()

