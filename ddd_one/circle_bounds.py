"""Exact maximum-extension check for DDD 3217. Coefficients are in Z[w], w of order 20."""
from collections import defaultdict
from itertools import combinations
from math import gcd

L = 5040
DEG = 8
ZERO = (0,) * DEG
ONE = (1,) + (0,) * (DEG - 1)
EXCEPTIONAL = {5: 1, 7: 3, 9: 2, 11: 4}


def times_w(a):
    # w^8 = w^6 - w^4 + w^2 - 1.
    b = [0, *a[:-1]]
    for i, sign in ((0,-1), (2,1), (4,-1), (6,1)):
        b[i] += sign * a[-1]
    return tuple(b)


POWERS = [ONE]
for _ in range(19):
    POWERS.append(times_w(POWERS[-1]))
assert times_w(POWERS[-1]) == ONE


def multiply(a, b):
    result = [0] * DEG
    for i, coefficient in enumerate(a):
        if coefficient:
            for j, other in enumerate(b):
                if other:
                    factor = coefficient * other
                    for k, x in enumerate(POWERS[i+j]):
                        result[k] += factor*x
    return tuple(result)


def primitive(p):
    while p and p[-1] == ZERO:
        p.pop()
    content = gcd(*(x for coefficient in p for x in coefficient))
    return [tuple(x//content for x in coefficient) for coefficient in p] if content else []


def remainder(p, q):
    # Pseudo-division needs no inversions of cyclotomic coefficients.
    while len(p) >= len(q):
        degree = len(p) - len(q)
        leading = p[-1]
        next_p = [multiply(q[-1], v) for v in p]
        for i, coefficient in enumerate(q):
            correction = multiply(leading, coefficient)
            next_p[i+degree] = tuple(x-y for x,y in zip(next_p[i+degree], correction))
        p = primitive(next_p)
    return p


def polynomial_gcd(p, q):
    p, q = list(p), list(q)
    while q:
        p, q = q, remainder(p, q)
    return p


def quadratic_resultant(p, q):
    c,b,a = p
    f,e,d = q
    def difference(x,y):
        return tuple(v-w for v,w in zip(x,y))
    af_cd = difference(multiply(a,f),multiply(c,d))
    ae_bd = difference(multiply(a,e),multiply(b,d))
    bf_ce = difference(multiply(b,f),multiply(c,e))
    return difference(multiply(af_cd,af_cd),multiply(ae_bd,bf_ce))


def reduced(records):
    circles = [r for r in records if any(r[0])]
    assert len(circles) == 6 and len({c[0] for c in circles}) == 1

    def contains(c, x):
        return all(len({(x[1][j]-c[1][j], (x[2][j]-c[2][j])%L)
                        for j in range(11) if c[0][j] == g}) == 1 for g in (0,1))

    fixed = [r for r in records if not any(r[0])]
    points = [r for r in fixed if not any(contains(c, r) for c in circles)]
    assert len(fixed) == 242 and len(points) == 2
    return circles + points


def orthogonality_polynomial(x, y):
    cells = defaultdict(list)
    for j in range(11):
        group = x[0][j] or y[0][j]
        e = (x[2][j] - y[2][j]) % L
        assert e % 252 == 0
        cells[group].append((x[1][j] - y[1][j], e//252))

    sums = defaultdict(lambda: [0]*DEG)
    def add(a, e, sign):
        for k, c in enumerate(POWERS[e%20]):
            sums[a][k] += sign*c

    if len(cells) == 1:
        for a,e in cells[0]:
            add(a,e,1)
    else:
        for group, terms in cells.items():
            for a,e in terms:
                for b,f in terms:
                    add(a-b,e-f,1 if group == 0 else -1)
    support = [a for a, coefficient in sums.items() if any(coefficient)]
    if not support:
        return []
    return primitive([tuple(sums[a]) for a in range(min(support),max(support)+1)])


def check_sharp(data, old, fixed_pair, cover, quotient):
    free = defaultdict(list)
    specified = defaultdict(set)
    for line in (data/'extensions.txt').read_text().splitlines():
        row = tuple(map(int,line.split()))
        h,z = row[:2]
        if z == -1:
            free[h].append((row[2:13],row[13:24],row[24:]))
        else:
            specified[h].add(z)
    assert set(free) == set(range(16))-set(EXCEPTIONAL)

    for filename in ('graphs.txt','audit_graphs.txt'):
        graphs = [tuple(map(int,s.split())) for s in (data/filename).read_text().splitlines()]
        fixed = [r for r in graphs if r[1]>=0]
        assert len(fixed) == 496 and all(r[-1]<=2 for r in fixed)
        assert all(r[-1]==1 for r in fixed if r[0] in EXCEPTIONAL)

    witnesses = {}
    zero_free_triples = 0
    monomial_exclusions = 0
    for h, records in sorted(free.items()):
        candidates = reduced(records)
        pair_polys = {(i,j): orthogonality_polynomial(candidates[i], candidates[j])
                      for i,j in combinations(range(8),2)}
        for i,j,k in combinations(range(8),3):
            p,q,r = pair_polys[i,j],pair_polys[i,k],pair_polys[j,k]
            if any(len(v)==1 for v in (p,q,r)):
                monomial_exclusions += 1
            else:
                quadratics = [v for v in (p,q,r) if v]
                assert len(quadratics)==2 and all(len(v)==3 for v in quadratics)
                assert quadratic_resultant(*quadratics)!=ZERO
            result = polynomial_gcd(polynomial_gcd(p,q),r)
            assert len(result)==1 and result[0]!=ZERO
            zero_free_triples += 1

        # Set the circle phases equal to one. Exhibit and test the five
        # complete rows as literal monomials, with z still free.
        plain = [tuple((a,0,e) for a,e in zip(row[1],row[2])) for row in candidates]
        witness = next(( (i,j) for i,j in combinations(range(8),2)
                        if cover(quotient(plain[i],plain[j])) ),None)
        assert witness is not None
        i,j = witness
        five = [((0,0,0),)*11, fixed_pair, old[h], plain[i], plain[j]]
        assert all(cover(quotient(x,y)) for x,y in combinations(five,2))
        witnesses[h] = {'circle_or_point_indices':[i,j],
                        'added_rows':[[[a,e//252] for a,g,e in row] for row in five[-2:]]}
        assert all(e%252==0 for row in five[-2:] for a,g,e in row)

    assert zero_free_triples == 672 and monomial_exclusions == 624
    coordinates = {1:(11,3,19,15,15,7),2:(15,3,11,19,15,7),
                   3:(11,7,19,15,3,15),4:(15,7,11,19,3,15)}
    for row in old:
        assert [j for j,(d,g,e) in enumerate(row) if d] == [6,7]
        assert tuple(row[j][2]//252 for j in (0,1,6,7,9)) == (0,5,0,10,10)
    exceptional = {}
    for h,r in EXCEPTIONAL.items():
        assert tuple(old[h][j][2]//252 for j in (2,3,4,5,8,10)) == coordinates[r]
        values = [m for m in range(60) if m%5!=0 and (m%5!=r or m%3==0)]
        assert len(values)==40 and specified[h]=={84*m for m in values}
        exceptional[h] = {'r':r,'exponents_mod_60':values}
    return {'uniform_five_row_families':witnesses,'exceptional_only_families':exceptional,
            'circle_or_point_triples':zero_free_triples,'monomial_exclusions':monomial_exclusions,
            'other_constant_gcds':zero_free_triples-monomial_exclusions,
            'maximum_total_rows':5}
