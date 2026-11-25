"""Exact value-multiset overlap audit for N=11 cycle profiles.

A p-cycle is the value set lambda*mu_p, with each value counted once. Two
cycles of different prime lengths share at most one value; two p-cycles that
meet have the same value set. For an overlap identity, relative phases along
each connected cycle-intersection component are products of roots whose
orders divide the LCM M. After normalizing one phase per component, an
exhaustive finite search over M-th-root phase classes decides existence.
"""
from collections import Counter, defaultdict
from itertools import combinations_with_replacement
from math import lcm

PROFILES = {
    "A": (11,),
    "B": (7, 2, 2),
    "C": (5, 3, 3),
    "D": (5, 2, 2, 2),
    "E": (3, 3, 3, 2),
    "F": (3, 2, 2, 2, 2),
}

def lcm_profiles(a, b):
    return lcm(*PROFILES[a], *PROFILES[b])

def cycle_options(p, M):
    """Distinct rotated p-cycles in the M-th roots, as multiplicity vectors."""
    out = []
    for a in range(M // p):
        v = [0] * M
        for k in range(p):
            v[(a + k * (M // p)) % M] += 1
        out.append(tuple(v))
    return tuple(out)

def add_vectors(vectors, M):
    total = [0] * M
    for v in vectors:
        for i, n in enumerate(v):
            total[i] += n
    return tuple(total)

def profile_multisets(profile, M):
    """Map each value-multiplicity vector to one cycle-phase witness."""
    out = {}
    ps = sorted(set(profile))
    def rec(k, chosen):
        if k == len(ps):
            out.setdefault(add_vectors(chosen, M), tuple(chosen))
            return
        p = ps[k]
        options = cycle_options(p, M)
        for pick in combinations_with_replacement(options, profile.count(p)):
            rec(k + 1, chosen + list(pick))
    rec(0, [])
    return out

def compact_witness(profile, M, witness):
    idx = 0
    result = {}
    for p in sorted(set(profile)):
        phases = []
        options = cycle_options(p, M)
        for _ in range(profile.count(p)):
            phases.append(options.index(witness[idx]))
            idx += 1
        result[p] = tuple(phases)
    return result

def vector_from_phases(phases, M):
    return add_vectors([cycle_options(p, M)[a]
                        for p, values in phases.items() for a in values], M)

def remove_common_cycles(phases_a, phases_b):
    common = {p: Counter(phases_a.get(p, ())) & Counter(phases_b.get(p, ()))
              for p in set(phases_a) | set(phases_b)}
    common = {p: c for p, c in common.items() if c}
    def residual(phases):
        out = {}
        for p, values in phases.items():
            left = Counter(values)
            left.subtract(common.get(p, Counter()))
            vals = tuple(sorted(x for x, n in left.items() for _ in range(max(0, n))))
            if vals:
                out[p] = vals
        return out
    return common, residual(phases_a), residual(phases_b)

# A cycle of prime length p, if not identical to an opposite-side cycle,
# must meet at least p distinct opposite-side cycles. This leaves exactly
# these four distinct-profile pairs for exact phase enumeration.
CANDIDATES = [("C", "D"), ("C", "F"), ("D", "F"), ("E", "F")]

if __name__ == "__main__":
    print("EXACT DISTINCT-PROFILE OVERLAP AUDIT")
    print("Candidate pairs after the cycle-count bound:", ", ".join(a+b for a,b in CANDIDATES))
    for a, b in CANDIDATES:
        M = lcm_profiles(a, b)
        left = profile_multisets(PROFILES[a], M)
        right = profile_multisets(PROFILES[b], M)
        common_vectors = left.keys() & right.keys()
        print(f"\n{a} ({'+'.join(map(str, PROFILES[a]))}) vs "
              f"{b} ({'+'.join(map(str, PROFILES[b]))}), M={M}: "
              f"{len(left)} vs {len(right)} distinct multiset vectors; "
              f"common={len(common_vectors)}")
        pattern_ok = True
        sample = None
        for value_vector in common_vectors:
            pa = compact_witness(PROFILES[a], M, left[value_vector])
            pb = compact_witness(PROFILES[b], M, right[value_vector])
            shared, ra, rb = remove_common_cycles(pa, pb)
            residual_a = vector_from_phases(ra, M)
            residual_b = vector_from_phases(rb, M)
            six_cycles = cycle_options(6, M) if M % 6 == 0 else ()
            if (residual_a != residual_b or residual_a not in six_cycles):
                pattern_ok = False
            if sample is None:
                sample = (pa, pb, shared)
        if common_vectors:
            print("  every common vector = shared prime cycles + one rotated mu_6:", pattern_ok)
            print("  sample phase classes (phase a means omega_M^a * mu_p):")
            print("   ", a, sample[0])
            print("   ", b, sample[1])
            print("  shared cycles:", sample[2])
        else:
            print("  no overlap at any phase assignment")
