"""Exact computation for the all-11-cycle stratum of regular 3x11 matrices.

After dephasing the first row, take
    H_pi = [ 1,...,1 ; zeta^j ; zeta^(-pi(j)) ],  j in Z/11Z,
where zeta is a primitive 11th root of unity. All three row quotients are
single rotated 11-cycles exactly when pi and j+pi(j) are permutations.
The script counts these complete mappings with dynamic programming and checks
the count independently by exhaustive permutation enumeration.
"""
from functools import lru_cache
from itertools import permutations

P = 11

def count_by_dynamic_programming(p=P):
    """Count complete mappings pi with pi(0)=0."""
    @lru_cache(None)
    def count(j, used_values, used_sums):
        if j == p:
            return 1
        if j == 0:
            return count(1, 1, 1)  # pi(0)=0 and 0+pi(0)=0
        total = 0
        for y in range(1, p):
            if (used_values >> y) & 1:
                continue
            s = (j + y) % p
            if (used_sums >> s) & 1:
                continue
            total += count(j + 1, used_values | (1 << y),
                           used_sums | (1 << s))
        return total
    return count(0, 0, 0)

def count_by_permutation_enumeration(p=P):
    """Independent exhaustive check over the (p-1)! normalized maps."""
    total = 0
    for tail in permutations(range(1, p)):
        sums = {0}
        for j, y in enumerate(tail, start=1):
            s = (j + y) % p
            if s in sums:
                break
            sums.add(s)
        else:
            total += 1
    return total

def main():
    dp = count_by_dynamic_programming()
    brute = count_by_permutation_enumeration()
    print(f"normalized complete mappings pi(0)=0: {dp}")
    print(f"independent exhaustive count: {brute}")
    print(f"counts agree: {dp == brute}")
    print(f"complete mappings before normalizing the v-row phase: {P * dp}")
    example = (0, 1, 2, 3, 5, 7, 4, 9, 6, 10, 8)
    sums = tuple((j + example[j]) % P for j in range(P))
    print(f"example pi: {example}")
    print(f"j+pi(j) mod 11: {sums}")

if __name__ == "__main__":
    main()
