# Code made for SageMath version 9.3
from math import ceil
from datetime import datetime
import time
# Setting up the symbolic variables
P.<r, s> = QQ[]
rs = [r, s]
r_one = rs.copy()
r_one.insert(0, 1)

U.<f0, f1, f2, f3, f4, f5> = QQ[]
f = U.gens()


# Returns f_n(x)
def n_angle(n, x):
    m = n // 2
    result = x^n
    for k in range(1, m + 1):
        p = binomial(n, 2 * k) * (x*x - 1)^k * x^(n - 2*k)
        result += p
    return result


# p is a polynomial or an integer
def find_unit(p):
    if isinstance(p, sage.rings.integer.Integer) or p.is_constant():
        unit_p = p
    else:
        unit_p = p.factor().unit()
    if isinstance(unit_p, sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular):
        unit_p = unit_p.constant_coefficient()
    return unit_p


# p is a fraction of 2 polynomials (they may be constant)
def reduce_fraction(p):
    num = p.numerator()
    unit_num = find_unit(num)
    den = p.denominator()
    unit_den = find_unit(den)
    unit_gcd = unit_num.gcd(unit_den)
    new_num = num / unit_gcd
    new_den = den / unit_gcd
    if unit_den < 0:
        return -new_num / (-new_den)
    else:
        return new_num / new_den


# Returns \cos(\widehat{s_1s_2s_3})
def get_angle(s1, s2, s3):
    return ((s1 + s3)^2 - (s1 + s2)^2 - (s2 + s3)^2) / (-2 * (s1 + s2) * (s2 + s3))


# Returns all the angles in the around-corona list_radii
def get_angles(list_radii, around):
    n = len(list_radii)
    result = {}
    for i in range(len(list_radii)):
        a = [list_radii[i], list_radii[(i + 1) % n]]
        a = sorted(a, key=lambda j: str(j))
        a.insert(1, around)
        lexi_angle = tuple(a)
        if lexi_angle in result:
            result[lexi_angle] += 1
        else:
            result[lexi_angle] = 1
    return result


# Precomputed formulas for Definition 3.9
addFormulas = [0, f0 - 1, f0 - f1,
               -2*f0*f1*f2 + f0^2 + f1^2 + f2^2 - 1,
               4*f0^2*f1^2*f2^2 - 4*f0^3*f1*f2*f3 - 4*f0*f1^3*f2*f3 - 4*f0*f1*f2^3*f3 + 4*f0^2*f1^2*f3^2 +
               4*f0^2*f2^2*f3^2 + 4*f1^2*f2^2*f3^2 - 4*f0*f1*f2*f3^3 + f0^4 - 2*f0^2*f1^2 + f1^4 - 2*f0^2*f2^2 -
               2*f1^2*f2^2 + f2^4 + 8*f0*f1*f2*f3 - 2*f0^2*f3^2 - 2*f1^2*f3^2 - 2*f2^2*f3^2 + f3^4,
               -128*f0^3*f1^3*f2^3*f3^3*f4^3 + 64*f0^4*f1^4*f2^2*f3^2*f4^2 + 64*f0^4*f1^2*f2^4*f3^2*f4^2 +
               64*f0^2*f1^4*f2^4*f3^2*f4^2 + 64*f0^4*f1^2*f2^2*f3^4*f4^2 + 64*f0^2*f1^4*f2^2*f3^4*f4^2 +
               64*f0^2*f1^2*f2^4*f3^4*f4^2 + 64*f0^4*f1^2*f2^2*f3^2*f4^4 + 64*f0^2*f1^4*f2^2*f3^2*f4^4 +
               64*f0^2*f1^2*f2^4*f3^2*f4^4 + 64*f0^2*f1^2*f2^2*f3^4*f4^4 - 32*f0^5*f1^3*f2^3*f3*f4 -
               32*f0^3*f1^5*f2^3*f3*f4 - 32*f0^3*f1^3*f2^5*f3*f4 - 32*f0^5*f1^3*f2*f3^3*f4 - 32*f0^3*f1^5*f2*f3^3*f4 -
               32*f0^5*f1*f2^3*f3^3*f4 - 32*f0*f1^5*f2^3*f3^3*f4 - 32*f0^3*f1*f2^5*f3^3*f4 - 32*f0*f1^3*f2^5*f3^3*f4 -
               32*f0^3*f1^3*f2*f3^5*f4 - 32*f0^3*f1*f2^3*f3^5*f4 - 32*f0*f1^3*f2^3*f3^5*f4 - 32*f0^5*f1^3*f2*f3*f4^3 -
               32*f0^3*f1^5*f2*f3*f4^3 - 32*f0^5*f1*f2^3*f3*f4^3 - 32*f0*f1^5*f2^3*f3*f4^3 - 32*f0^3*f1*f2^5*f3*f4^3 -
               32*f0*f1^3*f2^5*f3*f4^3 - 32*f0^5*f1*f2*f3^3*f4^3 - 32*f0*f1^5*f2*f3^3*f4^3 - 32*f0*f1*f2^5*f3^3*f4^3 -
               32*f0^3*f1*f2*f3^5*f4^3 - 32*f0*f1^3*f2*f3^5*f4^3 - 32*f0*f1*f2^3*f3^5*f4^3 - 32*f0^3*f1^3*f2*f3*f4^5 -
               32*f0^3*f1*f2^3*f3*f4^5 - 32*f0*f1^3*f2^3*f3*f4^5 - 32*f0^3*f1*f2*f3^3*f4^5 - 32*f0*f1^3*f2*f3^3*f4^5 -
               32*f0*f1*f2^3*f3^3*f4^5 + 16*f0^4*f1^4*f2^4 + 16*f0^6*f1^2*f2^2*f3^2 + 16*f0^2*f1^6*f2^2*f3^2 +
               16*f0^2*f1^2*f2^6*f3^2 + 16*f0^4*f1^4*f3^4 + 16*f0^4*f2^4*f3^4 + 16*f1^4*f2^4*f3^4 +
               16*f0^2*f1^2*f2^2*f3^6 + 16*f0^6*f1^2*f2^2*f4^2 + 16*f0^2*f1^6*f2^2*f4^2 + 16*f0^2*f1^2*f2^6*f4^2
               + 16*f0^6*f1^2*f3^2*f4^2 + 16*f0^2*f1^6*f3^2*f4^2 + 16*f0^6*f2^2*f3^2*f4^2 - 144*f0^4*f1^2*f2^2*f3^2*f4^2
               - 144*f0^2*f1^4*f2^2*f3^2*f4^2 + 16*f1^6*f2^2*f3^2*f4^2 - 144*f0^2*f1^2*f2^4*f3^2*f4^2
               + 16*f0^2*f2^6*f3^2*f4^2 + 16*f1^2*f2^6*f3^2*f4^2 - 144*f0^2*f1^2*f2^2*f3^4*f4^2 + 16*f0^2*f1^2*f3^6*f4^2
               + 16*f0^2*f2^2*f3^6*f4^2 + 16*f1^2*f2^2*f3^6*f4^2 + 16*f0^4*f1^4*f4^4 + 16*f0^4*f2^4*f4^4
               + 16*f1^4*f2^4*f4^4 - 144*f0^2*f1^2*f2^2*f3^2*f4^4 + 16*f0^4*f3^4*f4^4 + 16*f1^4*f3^4*f4^4
               + 16*f2^4*f3^4*f4^4 + 16*f0^2*f1^2*f2^2*f4^6 + 16*f0^2*f1^2*f3^2*f4^6 + 16*f0^2*f2^2*f3^2*f4^6
               + 16*f1^2*f2^2*f3^2*f4^6 - 8*f0^7*f1*f2*f3*f4 + 40*f0^5*f1^3*f2*f3*f4 + 40*f0^3*f1^5*f2*f3*f4
               - 8*f0*f1^7*f2*f3*f4 + 40*f0^5*f1*f2^3*f3*f4 + 112*f0^3*f1^3*f2^3*f3*f4 + 40*f0*f1^5*f2^3*f3*f4
               + 40*f0^3*f1*f2^5*f3*f4 + 40*f0*f1^3*f2^5*f3*f4 - 8*f0*f1*f2^7*f3*f4 + 40*f0^5*f1*f2*f3^3*f4
               + 112*f0^3*f1^3*f2*f3^3*f4 + 40*f0*f1^5*f2*f3^3*f4 + 112*f0^3*f1*f2^3*f3^3*f4 + 112*f0*f1^3*f2^3*f3^3*f4
               + 40*f0*f1*f2^5*f3^3*f4 + 40*f0^3*f1*f2*f3^5*f4 + 40*f0*f1^3*f2*f3^5*f4 + 40*f0*f1*f2^3*f3^5*f4
               - 8*f0*f1*f2*f3^7*f4 + 40*f0^5*f1*f2*f3*f4^3 + 112*f0^3*f1^3*f2*f3*f4^3 + 40*f0*f1^5*f2*f3*f4^3
               + 112*f0^3*f1*f2^3*f3*f4^3 + 112*f0*f1^3*f2^3*f3*f4^3 + 40*f0*f1*f2^5*f3*f4^3 + 112*f0^3*f1*f2*f3^3*f4^3
               + 112*f0*f1^3*f2*f3^3*f4^3 + 112*f0*f1*f2^3*f3^3*f4^3 + 40*f0*f1*f2*f3^5*f4^3 + 40*f0^3*f1*f2*f3*f4^5
               + 40*f0*f1^3*f2*f3*f4^5 + 40*f0*f1*f2^3*f3*f4^5 + 40*f0*f1*f2*f3^3*f4^5 - 8*f0*f1*f2*f3*f4^7
               - 8*f0^6*f1^2*f2^2 - 16*f0^4*f1^4*f2^2 - 8*f0^2*f1^6*f2^2 - 16*f0^4*f1^2*f2^4 - 16*f0^2*f1^4*f2^4
               - 8*f0^2*f1^2*f2^6 - 8*f0^6*f1^2*f3^2 - 16*f0^4*f1^4*f3^2 - 8*f0^2*f1^6*f3^2 - 8*f0^6*f2^2*f3^2
               - 24*f0^4*f1^2*f2^2*f3^2 - 24*f0^2*f1^4*f2^2*f3^2 - 8*f1^6*f2^2*f3^2 - 16*f0^4*f2^4*f3^2
               - 24*f0^2*f1^2*f2^4*f3^2 - 16*f1^4*f2^4*f3^2 - 8*f0^2*f2^6*f3^2 - 8*f1^2*f2^6*f3^2 - 16*f0^4*f1^2*f3^4
               - 16*f0^2*f1^4*f3^4 - 16*f0^4*f2^2*f3^4 - 24*f0^2*f1^2*f2^2*f3^4 - 16*f1^4*f2^2*f3^4 - 16*f0^2*f2^4*f3^4
               - 16*f1^2*f2^4*f3^4 - 8*f0^2*f1^2*f3^6 - 8*f0^2*f2^2*f3^6 - 8*f1^2*f2^2*f3^6 - 8*f0^6*f1^2*f4^2
               - 16*f0^4*f1^4*f4^2 - 8*f0^2*f1^6*f4^2 - 8*f0^6*f2^2*f4^2 - 24*f0^4*f1^2*f2^2*f4^2
               - 24*f0^2*f1^4*f2^2*f4^2 - 8*f1^6*f2^2*f4^2 - 16*f0^4*f2^4*f4^2 - 24*f0^2*f1^2*f2^4*f4^2
               - 16*f1^4*f2^4*f4^2 - 8*f0^2*f2^6*f4^2 - 8*f1^2*f2^6*f4^2 - 8*f0^6*f3^2*f4^2 - 24*f0^4*f1^2*f3^2*f4^2
               - 24*f0^2*f1^4*f3^2*f4^2 - 8*f1^6*f3^2*f4^2 - 24*f0^4*f2^2*f3^2*f4^2 + 192*f0^2*f1^2*f2^2*f3^2*f4^2
               - 24*f1^4*f2^2*f3^2*f4^2 - 24*f0^2*f2^4*f3^2*f4^2 - 24*f1^2*f2^4*f3^2*f4^2 - 8*f2^6*f3^2*f4^2
               - 16*f0^4*f3^4*f4^2 - 24*f0^2*f1^2*f3^4*f4^2 - 16*f1^4*f3^4*f4^2 - 24*f0^2*f2^2*f3^4*f4^2
               - 24*f1^2*f2^2*f3^4*f4^2 - 16*f2^4*f3^4*f4^2 - 8*f0^2*f3^6*f4^2 - 8*f1^2*f3^6*f4^2 - 8*f2^2*f3^6*f4^2
               - 16*f0^4*f1^2*f4^4 - 16*f0^2*f1^4*f4^4 - 16*f0^4*f2^2*f4^4 - 24*f0^2*f1^2*f2^2*f4^4 - 16*f1^4*f2^2*f4^4
               - 16*f0^2*f2^4*f4^4 - 16*f1^2*f2^4*f4^4 - 16*f0^4*f3^2*f4^4 - 24*f0^2*f1^2*f3^2*f4^4 - 16*f1^4*f3^2*f4^4
               - 24*f0^2*f2^2*f3^2*f4^4 - 24*f1^2*f2^2*f3^2*f4^4 - 16*f2^4*f3^2*f4^4 - 16*f0^2*f3^4*f4^4
               - 16*f1^2*f3^4*f4^4 - 16*f2^2*f3^4*f4^4 - 8*f0^2*f1^2*f4^6 - 8*f0^2*f2^2*f4^6 - 8*f1^2*f2^2*f4^6
               - 8*f0^2*f3^2*f4^6 - 8*f1^2*f3^2*f4^6 - 8*f2^2*f3^2*f4^6 - 24*f0^5*f1*f2*f3*f4 - 112*f0^3*f1^3*f2*f3*f4
               - 24*f0*f1^5*f2*f3*f4 - 112*f0^3*f1*f2^3*f3*f4 - 112*f0*f1^3*f2^3*f3*f4 - 24*f0*f1*f2^5*f3*f4
               - 112*f0^3*f1*f2*f3^3*f4 - 112*f0*f1^3*f2*f3^3*f4 - 112*f0*f1*f2^3*f3^3*f4 - 24*f0*f1*f2*f3^5*f4
               - 112*f0^3*f1*f2*f3*f4^3 - 112*f0*f1^3*f2*f3*f4^3 - 112*f0*f1*f2^3*f3*f4^3 - 112*f0*f1*f2*f3^3*f4^3
               - 24*f0*f1*f2*f3*f4^5 + f0^8 + 4*f0^6*f1^2 + 6*f0^4*f1^4 + 4*f0^2*f1^6 + f1^8 + 4*f0^6*f2^2
               + 28*f0^4*f1^2*f2^2 + 28*f0^2*f1^4*f2^2 + 4*f1^6*f2^2 + 6*f0^4*f2^4 + 28*f0^2*f1^2*f2^4 + 6*f1^4*f2^4
               + 4*f0^2*f2^6 + 4*f1^2*f2^6 + f2^8 + 4*f0^6*f3^2 + 28*f0^4*f1^2*f3^2 + 28*f0^2*f1^4*f3^2 + 4*f1^6*f3^2
               + 28*f0^4*f2^2*f3^2 + 40*f0^2*f1^2*f2^2*f3^2 + 28*f1^4*f2^2*f3^2 + 28*f0^2*f2^4*f3^2 + 28*f1^2*f2^4*f3^2
               + 4*f2^6*f3^2 + 6*f0^4*f3^4 + 28*f0^2*f1^2*f3^4 + 6*f1^4*f3^4 + 28*f0^2*f2^2*f3^4 + 28*f1^2*f2^2*f3^4
               + 6*f2^4*f3^4 + 4*f0^2*f3^6 + 4*f1^2*f3^6 + 4*f2^2*f3^6 + f3^8 + 4*f0^6*f4^2 + 28*f0^4*f1^2*f4^2
               + 28*f0^2*f1^4*f4^2 + 4*f1^6*f4^2 + 28*f0^4*f2^2*f4^2 + 40*f0^2*f1^2*f2^2*f4^2 + 28*f1^4*f2^2*f4^2
               + 28*f0^2*f2^4*f4^2 + 28*f1^2*f2^4*f4^2 + 4*f2^6*f4^2 + 28*f0^4*f3^2*f4^2 + 40*f0^2*f1^2*f3^2*f4^2
               + 28*f1^4*f3^2*f4^2 + 40*f0^2*f2^2*f3^2*f4^2 + 40*f1^2*f2^2*f3^2*f4^2 + 28*f2^4*f3^2*f4^2
               + 28*f0^2*f3^4*f4^2 + 28*f1^2*f3^4*f4^2 + 28*f2^2*f3^4*f4^2 + 4*f3^6*f4^2 + 6*f0^4*f4^4 + 28*f0^2*f1^2*f4^4
               + 6*f1^4*f4^4 + 28*f0^2*f2^2*f4^4 + 28*f1^2*f2^2*f4^4 + 6*f2^4*f4^4 + 28*f0^2*f3^2*f4^4 + 28*f1^2*f3^2*f4^4
               + 28*f2^2*f3^2*f4^4 + 6*f3^4*f4^4 + 4*f0^2*f4^6 + 4*f1^2*f4^6 + 4*f2^2*f4^6 + 4*f3^2*f4^6 + f4^8
               + 72*f0^3*f1*f2*f3*f4 + 72*f0*f1^3*f2*f3*f4 + 72*f0*f1*f2^3*f3*f4 + 72*f0*f1*f2*f3^3*f4
               + 72*f0*f1*f2*f3*f4^3 - 4*f0^6 - 12*f0^4*f1^2 - 12*f0^2*f1^4 - 4*f1^6 - 12*f0^4*f2^2
               - 32*f0^2*f1^2*f2^2 - 12*f1^4*f2^2 - 12*f0^2*f2^4 - 12*f1^2*f2^4 - 4*f2^6 - 12*f0^4*f3^2
               - 32*f0^2*f1^2*f3^2 - 12*f1^4*f3^2 - 32*f0^2*f2^2*f3^2 - 32*f1^2*f2^2*f3^2 - 12*f2^4*f3^2 - 12*f0^2*f3^4
               - 12*f1^2*f3^4 - 12*f2^2*f3^4 - 4*f3^6 - 12*f0^4*f4^2 - 32*f0^2*f1^2*f4^2 - 12*f1^4*f4^2
               - 32*f0^2*f2^2*f4^2 - 32*f1^2*f2^2*f4^2 - 12*f2^4*f4^2 - 32*f0^2*f3^2*f4^2 - 32*f1^2*f3^2*f4^2
               - 32*f2^2*f3^2*f4^2 - 12*f3^4*f4^2 - 12*f0^2*f4^4 - 12*f1^2*f4^4 - 12*f2^2*f4^4 - 12*f3^2*f4^4
               - 4*f4^6 - 40*f0*f1*f2*f3*f4 + 6*f0^4 + 12*f0^2*f1^2 + 6*f1^4 + 12*f0^2*f2^2 + 12*f1^2*f2^2 + 6*f2^4
               + 12*f0^2*f3^2 + 12*f1^2*f3^2 + 12*f2^2*f3^2 + 6*f3^4 + 12*f0^2*f4^2 + 12*f1^2*f4^2 + 12*f2^2*f4^2
               + 12*f3^2*f4^2 + 6*f4^4 - 4*f0^2 - 4*f1^2 - 4*f2^2 - 4*f3^2 - 4*f4^2 + 1]


# Returns the polynomial for each angle, which is f_n(\cos(\widehat{angle})) where n is the amount of times
# the angle appears in the r-corona
def get_angle_formulae(angles):
    angle_formulae = []
    for angle in angles:
        formula = n_angle(angles[angle], get_angle(*angle))
        angle_formulae.append(formula)
    return angle_formulae


# Moves each element of cor i places to the left
def permutate_cyclically(cor, i):
    new_cor = []
    for j in range(len(cor)):
        new_cor.append(cor[(i + j) % len(cor)])
    return new_cor


sCoronasFernique = []
rCoronasFernique = []


# Checks if there exists a corona for a given list of angle multiplicities (the vector \vec{h})
def check_angle_coding(ks):
    return (ks[1] + ks[4]) % 2 == 0 and (ks[2] + ks[4]) % 2 == 0 and \
        (ks[0] == 0 or (ks[1] != 0 or ks[2] != 0) or ks[1] == ks[2] == ks[3] == ks[4] == ks[5] == 0) and \
        (ks[3] == 0 or (ks[1] != 0 or ks[4] != 0) or ks[1] == ks[2] == ks[0] == ks[4] == ks[5] == 0) and \
        (ks[5] == 0 or (ks[4] != 0 or ks[2] != 0) or ks[1] == ks[2] == ks[3] == ks[4] == ks[0] == 0)


# Gets all the coronas where needed is the amount of 1's, r's and s's we want
def get_coronas_recursive(needed, cur_list, permutations, current):
    if len(cur_list) == sum(needed):
        cyclic_permutations = []
        for i in range(len(cur_list)):
            cyclic_permutations.append(permutate_cyclically(cur_list, i))
            cyclic_permutations.append([])
            for j in range(len(cur_list)):
                cyclic_permutations[-1].append(cur_list[(i - j) % len(cur_list)])
        cyclic_permutations.sort()
        if cyclic_permutations[0] == cur_list:
            permutations.append(cur_list)
    else:
        for i in range(len(needed)):
            if current[i] < needed[i]:
                new_list = cur_list.copy()
                new_list.append(i)
                new_current = current.copy()
                new_current[i] += 1
                get_coronas_recursive(needed, new_list, permutations, new_current)
    return permutations


# Gets the index for \vec{h}
def get_index(a, b):
    if a == 0 and b == 0:
        return 0
    elif (a == 1 and b == 0) or (a == 0 and b == 1):
        return 1
    elif (a == 2 and b == 0) or (a == 0 and b == 2):
        return 2
    elif a == 1 and b == 1:
        return 3
    elif (a == 1 and b == 2) or (a == 2 and b == 1):
        return 4
    else:
        return 5


# Gets the multiplicity vector \vec{h} for a permutation
def get_coding(perm):
    cur_k = [0] * 6
    for i in range(len(perm)):
        cur_k[get_index(perm[i], perm[(i+1)%len(perm)])] += 1
    return cur_k


def get_coding2(perm):
    cur_k = [0] * 6
    for i in range(len(perm)):
        cur_k[get_index(r_one.index(perm[i]), r_one.index(perm[(i+1)%len(perm)]))] += 1
    return cur_k


# Gets all the coronas from a given multiplicity vector ks
def get_coronas_from_coding(ks):
    ones = ks[0] + (ks[1] + ks[2]) // 2
    rrs = (ks[1] + ks[4]) // 2 + ks[3]
    ss = (ks[2] + ks[4]) // 2 + ks[5]
    perms = get_coronas_recursive([ones, rrs, ss], [], [], [0, 0, 0])
    good_perms = []
    for perm in perms:
        cur_k = get_coding(perm)
        if cur_k == list(ks):
            good_perms.append([r_one[i] for i in perm])
    return good_perms


# Gets all the coronas for all multipicity vectors kss
def get_coronas_from_checked_coding(kss):
    output = []
    for ks in kss:
        output.extend(get_coronas_from_coding(ks))
    return output


def get_s_coronas_fernique():
    global sCoronasFernique
    sCoronasFernique = []
    s_angles = []
    # wrong_angles = []   # remove
    for k1 in range(6):
        for k2 in range(6 - k1):
            for k3 in range(6 - k1 - k2):
                for k4 in range(6 - k1 - k2 - k3):
                    for k5 in range(6 - k1 - k2 - k3 - k4):
                        for k6 in range(6 - k1 - k2 - k3 - k4 - k5):
                            if 6 < 3*k1 + 3*k2 + 3/2*k3 + 3*k4 + 3/2*k5 + k6 and \
                                    check_angle_coding((k1, k2, k3, k4, k5, k6)):
                                s_angles.append((k1, k2, k3, k4, k5, k6))
                            # else:
                            #     wrong_angles.append(''.join([str(_k) for _k in (k1, k2, k3, k4, k5, k6)]))
    # print('Wrong angles:', wrong_angles)
    sCoronasFernique = get_coronas_from_checked_coding(s_angles)
    s_coronas = [''.join([str(_k) for _k in s_corona]) for s_corona in sCoronasFernique]
    print(s_coronas)

get_s_coronas_fernique()


def make_permutations(length, cur_list, permutations, n):
    if len(cur_list) == length:
        cyclic_permutations = []
        for i in range(len(cur_list)):
            cyclic_permutations.append(permutate_cyclically(cur_list, i))
            cyclic_permutations.append([])
            for j in range(len(cur_list)):
                cyclic_permutations[-1].append(cur_list[(i - j) % len(cur_list)])
        cyclic_permutations.sort()
        if cyclic_permutations[0] == cur_list:
            permutations.append(cur_list)
    else:
        for i in range(n):
            new_list = cur_list.copy()
            new_list.append(i)
            make_permutations(length, new_list, permutations, n)
    return permutations


rAnglesFernique = []
margin = 10 ^ (-5)


def get_r_coronas_fernique(a):
    global rAnglesFernique
    rAnglesFernique = []
    u_a = arccos(1 / (1 + a))
    v_a = arccos(1 - (2*a*a / (1+a)^2))
    f_u = 3 / pi * u_a
    f_v = 3 / pi * v_a
    if_u = int(ceil(6 / f_u))
    if_v = int(ceil(6 / f_v))
    r_angles = []
    for k3 in range(if_u):
        for k5 in range(if_u - int(k3 * f_u)):
            for k6 in range(if_v - int(k3 * f_u + k5 * f_u)):
                for k1 in range(6 - int(k3 * f_u + k5 * f_u + k6 * f_v)):
                    for k2 in range(6 - int(k3 * f_u + k5 * f_u + k6 * f_v) - k1):
                        for k4 in range(6 - int(k3 * f_u + k5 * f_u + k6 * f_v) - k1 - k2):
                            if 6 < 3*k1 + 3/2*k2 + 3/2*k3 + k4 + k5 + k6 and \
                                    k1 + k2 + k4 + (k3+k5)/2 < 6 and check_angle_coding((k1, k2, k3, k4, k5, k6)):
                                r_angles.append((k1, k2, k3, k4, k5, k6))
    rAnglesFernique = r_angles


# Gets all the possible multiplicity vectors for a 1-corona, given a solution
def get_angles_around_one(solution, s_corona, fast=False):
    _r, _s = solution
    # get_angle(1, 1, 1) == 1/2
    cos_angles = [1/2, get_angle(1, 1, _r), get_angle(1, 1, _s),
                  get_angle(_r, 1, _r), get_angle(_r, 1, _s), get_angle(_s, 1, _s)]
    print(arccos(cos_angles[0]), type(arccos(cos_angles[0])))
    angles = [arccos(angle).numerical_approx() for angle in cos_angles]
    kss = []
    s_has_one = 1 in s_corona
    for k1 in range(int(2*pi / angles[0]) + 2):
        sub = k1*angles[0]
        for k2 in range(int((2*pi - sub) / angles[1]) + 2):
            sub = k1*angles[0] + k2*angles[1]
            for k3 in range(int((2*pi - sub) / angles[2]) + 2):
                if k3 > 0 and not s_has_one:
                    continue
                sub = k1*angles[0] + k2*angles[1] + k3*angles[2]
                for k4 in range(int((2*pi - sub) / angles[3]) + 2):
                    sub = k1*angles[0] + k2*angles[1] + k3*angles[2] + k4*angles[3]
                    for k5 in range(int((2*pi - sub) / angles[4]) + 2):
                        if k5 > 0 and not s_has_one:
                            continue
                        sub = k1*angles[0] + k2*angles[1] + k3*angles[2] + k4*angles[3] + k5*angles[4]
                        for k6 in range(int((2*pi - sub) / angles[5]) + 2):
                            if k6 > 0 and not s_has_one:
                                continue
                            sub = k1*angles[0] + k2*angles[1] + k3*angles[2] + k4*angles[3] + k5*angles[4] + k6*angles[5]
                            sub2 = sub - 2*pi
                            sub2 = sub2.numerical_approx()
                            if abs(sub2) < margin and check_angle_coding((k1, k2, k3, k4, k5, k6)):
                                kss.append((k1, k2, k3, k4, k5, k6))
                                if fast and (k1, k2, k3, k4, k5, k6) != (6, 0, 0, 0, 0, 0):
                                    return kss
    if kss[0] == (6, 0, 0, 0, 0, 0) and len(kss) == 1:
        return []
    return kss


# Gets all the 1-coronas for a solution
def get_coronas_around_one(solution):
    kss = get_angles_around_one(solution)
    return get_coronas_from_checked_coding(kss)


def verify_solutions(solutions, angles):
    verified_sols = []
    for sol, polys in solutions:
        cur_angle = 0
        for angle in angles:
            formula = get_angle(*angle)
            sub_formula = {rs[i]: sol[i] for i in range(len(sol))}
            cos_f = formula.subs(sub_formula)
            a = arccos(cos_f)
            cur_angle += a * angles[angle]
            print(type(cur_angle))
            print('Cur angle:', cur_angle)
        if bool(pi * 2 - margin < cur_angle) and bool(cur_angle < pi * 2 + margin):
            verified_sols.append((sol, polys))
    return verified_sols


# Returns the roots of a polynomial
def get_roots(poly):
    if poly == 0:
        return None
    poly = poly.factor()
    factors = list(poly)
    total_roots = []
    for factor, count in factors:
        if isinstance(factor, sage.rings.integer.Integer):
            continue
        if factor.nvariables() == 1 and len(factor.coefficients()) > 1:
            try:
                ufactor = factor.univariate_polynomial()
                roots = ufactor.roots(RIF)
                for root, mult in roots:
                    if 0 <= root <= 1 and root.upper() > 0 and root.lower() < 1:
                        total_roots.append(root)
            except Exception as e:
                roots = ufactor.roots(CIF)
                for root, mult in roots:
                    if abs(root.imag()).lower() == 0 and root.real().lower() >= 0 and root.real().upper() <= 1:
                        total_roots.append(root)
    return total_roots


def file_print(st, file):
    print(st)
    if file is not None:
        file.write(st)
        file.write('\n')


def get_factors_angles(angles):
    add_formula = addFormulas[len(angles)]
    angle_formulae = get_angle_formulae(angles)
    for i in range(len(angle_formulae)):
        angle_formulae[i] = reduce_fraction(angle_formulae[i])
    add_formula = add_formula.subs({f[i]: formula for i, formula in enumerate(angle_formulae)})
    num = add_formula.numerator()
    num = num.factor()
    factors = list(num)
    good_factors = []
    for factor, count in factors:
        if factor.nvariables() > 0 and len(factor.coefficients()) > 1:
            good_factors.append(factor)
    return good_factors


def get_possible_radii_angles(s_corona, r_angles, file=None, fast=False, extra_poly=None):
    solutions = []
    s_angles = get_angles(s_corona, s)
    if len(r_angles) >= 6 or len(s_angles) >= 6:
        return True
    good_factors_s = get_factors_angles(s_angles)
    good_factors_r = get_factors_angles(r_angles)
    for i in range(len(good_factors_s)):
        for j in range(len(good_factors_r)):
            p1 = good_factors_s[i]
            p2 = good_factors_r[j]
            res1 = p1.resultant(p2, r)
            res2 = p2.resultant(p1, s)
            if p1 == r - s or p1 == -r + s or p2 == r - s or p2 == -r + s:
                continue  # Since r > s, we can have no solutions where r = s.
            roots_s = get_roots(res1)
            roots_r = get_roots(res2)
            if roots_s is None or roots_r is None:
                if extra_poly is None:
                    possible = False
                    # If the resultant is 0, either polynomial is a multiple of the other.
                    if p1.divides(p2):
                        extra_poly = p2
                    else:
                        extra_poly = p1
                    for new_r_angles in get_2_radii_angles():
                        possible2 = get_possible_radii_angles(s_corona, new_r_angles, file, fast, extra_poly)
                        if possible2:
                            return True
                    return possible
                else:
                    return True
            for root_r in roots_r:
                for root_s in roots_s:
                    if root_s < root_r:
                        if extra_poly is None:
                            solutions.append(((root_r, root_s), (res1, res2)))
                        else:
                            x = extra_poly(root_r, root_s)
                            if -margin < x < margin:
                                solutions.append(((root_r, root_s), (res1, res2)))
    verified_sols = verify_solutions(solutions, s_angles)
    if len(verified_sols) > 0:
        verified_sols = verify_solutions(verified_sols, r_angles)
        for sol, polys in verified_sols:
            one_angles = get_angles_around_one(sol, s_corona, fast=fast)
            file_print('One-angles: {}'.format(one_angles), file)
            if len(one_angles) > 0:
                file_print('Radii: {}'.format(sol), file)
                return one_angles
    return []


sCoronasData = []


def match_s_pattern(l0, l1):
    if len(l0) != len(l1):
        return False
    for i in range(len(l0)):
        if (l0[i] == s and l1[i] != s) or (l0[i] != s and l1[i] == s):
            return False
    return True


def really_match_pattern(l0, l1):
    if len(l0) != len(l1):
        return False
    for i in range(len(l0)):
        if match_s_pattern(permutate_cyclically(l0, i), l1):
            return True
    return False


def create_s_data():
    global sCoronasData
    sCoronasData = [(0.701, [[1, 1, 1, 1, 1]]),
                    (0.637, [[1, 1, 1, 1, s]]),
                    (0.545, [[1, 1, 1, s, s]]),
                    (0.533, [[1, 1, s, 1, s]]),
                    (0.414, [[1, 1, 1, 1]]),
                    (0.386, [[1, 1, s, s, s]]),
                    (0.349, [[1, s, 1, s, s]]),
                    (0.280, [[1, 1, 1, s]]),
                    (0.154, [[1, 1, 1]]),
                    (0.101, [[1, 1, s, s]])]
    for s_c in sCoronasFernique:
        for i in range(len(sCoronasData)):
            if s_c in sCoronasData[i][1]:
                break
            elif really_match_pattern(s_c, sCoronasData[i][1][0]):
                sCoronasData[i][1].append(s_c)
                break


def get_alpha(s_corona):
    for i in range(len(sCoronasData)):
        if s_corona in sCoronasData[i][1]:
            return sCoronasData[i][0]
    return 0.101


def get_2_radii_angles():
    r_angless = []
    for i in range(len(sCoronasData)):
        r_angless.append(get_angles(list(map(lambda x: r if x == s else x, sCoronasData[i][1][0])), r))
    return r_angless


def get_total_pairs():
    total = 0
    for alpha, c_list in sCoronasData:
        get_r_coronas_fernique(alpha)
        total += len(c_list) * len(rAnglesFernique)
    return total


create_s_data()
get_total_pairs()


def get_possible_radii_coding(s_corona, r_coding, file=None, fast=False):
    if r in s_corona and r_coding[2] == 0 and r_coding[4] == 0 and r_coding[5] == 0:
        return False
    s_coding = get_coding2(s_corona)
    if (s_coding[1] > 0 and r_coding[2] == 0) or (s_coding[3] > 0 and r_coding[4] == 0) or \
       (s_coding[4] > 0 and r_coding[5] == 0) or (r_coding[2] > 0 and s_coding[1] == 0) or \
       (r_coding[4] > 0 and s_coding[3] == 0) or (r_coding[5] > 0 and s_coding[4] == 0):
        return False
    r_angles = {}
    for i in range(len(r_coding)):
        if r_coding[i] > 0:
            if i == 0:
                r_angles[(1, r, 1)] = r_coding[i]
            elif i == 1:
                r_angles[(1, r, r)] = r_coding[i]
            elif i == 2:
                r_angles[(1, r, s)] = r_coding[i]
            elif i == 3:
                r_angles[(r, r, r)] = r_coding[i]
            elif i == 4:
                r_angles[(r, r, s)] = r_coding[i]
            else:
                r_angles[(s, r, s)] = r_coding[i]
    r_coronas = get_coronas_from_coding(r_coding)
    if len(r_coronas) == 1:
        r_corona = r_coronas[0]
        for i in range(len(r_corona)):
            if r_corona[i] == s:
                sub_str = str(r_corona[(i-1) % len(r_corona)]) + 'r' + str(r_corona[(i+1) % len(r_corona)])
                s_corona_str = ''.join([str(j) for j in s_corona])
                if not(sub_str in s_corona_str * 2 or sub_str in reversed(s_corona_str * 2)):
                    return False
    file_print('Analysing s: {} with r_angles {}'.format(s_corona, r_coding), file)
    return get_possible_radii_angles(s_corona, r_angles, file=file, fast=fast)


def check_this_s_corona(s_corona, file=None, do_alpha=True):
    r_codings = []
    start = time.time()
    if do_alpha:
        get_r_coronas_fernique(get_alpha(s_corona))
    for i, r_angles in enumerate(rAnglesFernique):
        print('Analysing {}/{}...'.format(i+1, len(rAnglesFernique)))
        if get_possible_radii_coding(s_corona, r_angles, file=file, fast=True):
            r_codings.append(r_angles)
    end = time.time()
    print('Time:', end - start)
    print('R-codings:', r_codings)
    print('#R-codings:', len(r_codings))
    return r_codings


def get_latex_corona(s_corona):
    return ''.join([str(x).upper() for x in s_corona])


def do_everything():
    total_pairs = {}
    now = datetime.now()
    with open('output_{}_{}_{}.txt'.format(now.hour, now.minute, now.second), 'w') as file:
        alpha, s_coronas = sCoronasData[3]
        get_r_coronas_fernique(alpha)
        for j, s_corona in enumerate(s_coronas):
            print('---------------------------------')
            print('---------------------------------')
            print('ANALYSING S-CORONA {}/{}'.format(j + 1, len(s_coronas)))
            print('---------------------------------')
            print('---------------------------------')
            r_codings = check_this_s_corona(s_corona, file, False)
            total_pairs[get_latex_corona(s_corona)] = r_codings
        print('---------------------------------')
        print('---------------------------------')
        file_print('TOTAL PAIRS: {}'.format(total_pairs), file)
        file_print('#TOTAL PAIRS: {}'.format(len(total_pairs)), file)
        print('---------------------------------')
        print('---------------------------------')
    print('TOTAL PAIRS:', total_pairs)
    print('#TOTAL PAIRS:', len(total_pairs))


do_everything()