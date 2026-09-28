#!/usr/bin/env python3
"""Exact finite algebra checks accompanying III-W (Python 3, standard library).

Checks selected rational moment integrals, discrete overlap identities,
polynomial instances of the transverse quadratic identities, and the
algebra of the cosine-taper certificate. These are reproducibility checks,
NOT proofs of the arithmetic estimates, infinite-dimensional stability,
limit exchanges, spectral laws, or asymptotic theorems in the manuscript.
No files, network, external data, or third-party packages are used.
"""

from fractions import Fraction as F
from math import comb, factorial
from random import Random


# A polynomial is a list of rational coefficients in ascending order.
def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p or [F(0)]


def add(p, q):
    r = [F(0)] * max(len(p), len(q))
    for i, a in enumerate(p):
        r[i] += a
    for i, a in enumerate(q):
        r[i] += a
    return trim(r)


def scale(p, c):
    return trim([F(c) * a for a in p])


def sub(p, q):
    return add(p, scale(q, -1))


def mul(p, q):
    r = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i + j] += a * b
    return trim(r)


def power(p, n):
    r = [F(1)]
    for _ in range(n):
        r = mul(r, p)
    return r


def primitive(p):
    return [F(0)] + [a / (i + 1) for i, a in enumerate(p)]


def integral(p):
    return sum((a / (i + 1) for i, a in enumerate(p)), F(0))


def reflect(p):
    """p(1-x), exactly."""
    r = [F(0)]
    for i, a in enumerate(p):
        r = add(r, scale(power([F(1), F(-1)], i), a))
    return r


def norm_squared(p):
    return integral(mul(p, p))


def real_correlation(p):
    """Polynomial in r equal to integral_0^(1-r) p(x)p(x+r) dx."""
    rpoly = [F(0)]
    for i, a in enumerate(p):
        for j, b in enumerate(p):
            for k in range(j + 1):
                term = [F(0)] * (j - k) + power(
                    [F(1), F(-1)], i + k + 1
                )
                coefficient = a * b * F(comb(j, k), i + k + 1)
                rpoly = add(rpoly, scale(term, coefficient))
    return trim(rpoly)


def check_low_moments():
    # Exact polynomial integration after splitting into triangles.
    m2 = 2 * integral([F(0), F(1), F(-1)])
    imax = integral([F(0), F(0), F(0), F(1), F(-1)])
    # Integral over a,b>=0, a+b<=1 of a*b*(1-a-b).
    isum = F(factorial(1) * factorial(1) * factorial(1), factorial(5))
    assert (m2, imax, isum) == (F(1, 3), F(1, 20), F(1, 120))
    m4 = 4 * imax + 8 * isum
    xxyy = 2 * imax + 2 * isum
    xyxy = 4 * isum
    commutator = 2 * (xxyy - xyxy)
    assert (m4, xxyy, xyxy, commutator) == (
        F(4, 15), F(7, 60), F(1, 30), F(1, 6)
    )
    assert F(1, 2) - F(1, 12) == F(5, 12)
    assert F(1, 2) ** 2 * F(1, 12) == F(1, 48)
    print("PASS: moment integrals, mixed moments, commutator, base constants")


def check_discrete_identities():
    rng = Random(20260927)
    count = 0
    for d in range(1, 25):
        kernel = [F(r * (d - r)) for r in range(d)]
        assert sum(kernel, F(0)) / d == F(d * d - 1, 6)
        for r in range(d):
            second_difference = (
                kernel[(r + 1) % d] - 2 * kernel[r] + kernel[(r - 1) % d]
            )
            assert second_difference == -2 + (2 * d if r == 0 else 0)
        flat = sum(
            (j * F(d - j, d) ** 2 for j in range(1, d)), F(0)
        )
        assert flat == F(d * d - 1, 12)

        # Homogeneous physical-space form of the discrete deficit.
        for _ in range(3):
            v = [F(rng.randrange(-5, 6), rng.randrange(1, 5)) for _ in range(d)]
            density = [x * x for x in v]
            mass = sum(density, F(0))
            upper = F(d * d - 1, 12) * mass * mass
            objective = F(0)
            overlap = F(0)
            for j in range(1, d):
                products = [v[i] * v[i + j] for i in range(d - j)]
                objective += j * sum(products, F(0)) ** 2
                overlap += F(j, 2) * sum(
                    ((a - b) ** 2 for a in products for b in products), F(0)
                )
            kernel_energy = F(1, 2) * sum(
                (
                    kernel[(i - k) % d] * density[i] * density[k]
                    for i in range(d) for k in range(d)
                ),
                F(0),
            )
            assert upper - objective == upper - kernel_energy + overlap
            assert objective <= upper
            count += 1
    print(f"PASS: discrete kernel, flat equality, {count} exact deficit samples")


def direct_quadratic(p, imaginary=False):
    # Called only for mean-zero p.
    prim = primitive(p)
    if imaginary:
        linear = scale(add(prim, reflect(prim)), -1)
    else:
        linear = sub(reflect(prim), prim)
    correlation = real_correlation(p)
    return (
        norm_squared(p) / 6
        - integral(mul([F(0), F(1)], mul(linear, linear)))
        - 2 * integral(mul([F(0), F(1), F(-1)], correlation))
    )


def check_transverse_quadratic():
    rng = Random(618)
    count = 0
    for degree in range(1, 9):
        for _ in range(3):
            p = [F(rng.randrange(-4, 5), 3) for _ in range(degree + 1)]
            # Real tangent complement: integral a = 0.
            a = sub(p, [integral(p)])
            assert integral(a) == 0
            ap = primitive(a)
            ae = scale(add(ap, reflect(ap)), F(1, 2))
            ae_centered = sub(ae, [integral(ap)])
            assert direct_quadratic(a) == (
                norm_squared(a) / 6 + 2 * norm_squared(ae_centered)
            )
            # Imaginary tangent complement: integral b = integral x*b = 0.
            m0 = integral(p)
            m1 = integral(mul([F(0), F(1)], p))
            b = sub(p, [4 * m0 - 6 * m1, 12 * m1 - 6 * m0])
            assert integral(b) == 0
            assert integral(mul([F(0), F(1)], b)) == 0
            bp = primitive(b)
            bo = scale(sub(bp, reflect(bp)), F(1, 2))
            assert direct_quadratic(b, imaginary=True) == (
                norm_squared(b) / 6 + 2 * norm_squared(bo)
            )
            count += 2
    # A nonzero real polynomial direction attaining the quadratic gap.
    a = [F(1), F(-6), F(6)]
    assert integral(a) == 0
    assert direct_quadratic(a) == norm_squared(a) / 6
    assert F(1, 6) - F(7, 3) * F(1, 32) - F(1, 24 * 32**2) >= F(1, 12)
    print(f"PASS: {count} polynomial transverse identities and local-radius constant")


def check_taper_algebra():
    # y=n^2; R(y)=(y-1)^2/(y(y+1)).
    numerator = [F(1), F(-2), F(1)]
    denominator = [F(0), F(1), F(1)]
    numerator_derivative = [F(-2), F(2)]
    denominator_derivative = [F(1), F(2)]
    derivative_numerator = sub(
        mul(numerator_derivative, denominator),
        mul(numerator, denominator_derivative),
    )
    assert derivative_numerator == mul([F(-1), F(1)], [F(1), F(3)])
    assert F((4 - 1) ** 2, 4 * (4 + 1)) == F(9, 20)

    # z=pi^2: Khat(1)=(-24+epsilon*(4*z-3))/(48*z).
    assert sub([F(0), F(1, 12)], [F(1, 16)]) == scale(
        [F(-3), F(4)], F(1, 48)
    )
    assert F(-1, 2) == F(-24, 48)
    # Pressure is half the constant Fourier coefficient.
    k0 = [F(1, 6), F(-1, 2)]  # coefficients in epsilon/pi^2
    pressure = [F(1, 12), F(-1, 4)]
    assert scale(k0, F(1, 2)) == pressure
    print("PASS: cosine certificate endpoint algebra and pressure normalization")


def main():
    check_low_moments()
    check_discrete_identities()
    check_transverse_quadratic()
    check_taper_algebra()
    print("All finite checks passed. The manuscript supplies the analytic proofs.")


if __name__ == "__main__":
    main()
