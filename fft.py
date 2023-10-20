import numpy as np
import numpy.polynomial.polynomial as polynomial
import cmath
import math

def mult(poly1 : polynomial.Polynomial, poly2 : polynomial.Polynomial):
    res_deg = poly1.degree() + poly2.degree()
    n = next_power_of_2(res_deg + 1)
    poly1 = make_deg(poly1, n - 1)
    poly2 = make_deg(poly2, n - 1)
    omega = get_omega(n)
    points1 = fft(poly1, omega)
    points2 = fft(poly2, omega)
    points_res = [ y1 * y2 for (y1, y2) in zip(points1, points2)]
    poly_res = inv_fft(points_res, omega)
    return poly_res

def fft(poly : polynomial.Polynomial, w : complex):
    if (cmath.isclose(w, complex(1,0))):
        return [poly.coef[0]]
    n = poly.degree() + 1
    [even, odd] = split_poly(poly)
    s1 = fft(even, w * w)
    s2 = fft(odd, w * w)
    r = np.zeros(n, complex)
    for j in range(0, (n // 2)):
        r[j] = s1[j] + (w ** j) * s2[j]
        r[j + (n // 2)]  = s1[j] - (w ** j) * s2[j]
    return r

def inv_fft(points, w):
    poly = fft(polynomial.Polynomial(points), (1 / w))
    n = len(poly)
    res = [coef * (1 / n) for coef in poly]
    return polynomial.Polynomial(res)

def split_poly(poly : polynomial.Polynomial):
    assert(poly.degree() % 2 == 1)
    even = []
    odd = []
    for i in range(0, poly.degree() + 1, 2):
        even.append(poly.coef[i])
        odd.append(poly.coef[i + 1])
    return [polynomial.Polynomial(even), polynomial.Polynomial(odd)]

def make_deg(poly : polynomial.Polynomial, deg : int):
    if poly.degree() == deg:
        return poly.copy()
    if poly.degree() > deg:
        return poly.cutdeg(deg)
    else:
        coefs = np.concatenate((poly.coef, np.zeros(deg - poly.degree())))
        return polynomial.Polynomial(coefs)

def to_round_poly(poly : polynomial.Polynomial):
    real_poly = [np.real(c) for c in poly.coef]
    return polynomial.Polynomial(np.rint(real_poly))

def align_polynomials(poly1 : polynomial.Polynomial, poly2 : polynomial.Polynomial):
    req_deg = max(poly1.degree(), poly2.degree())
    poly1 = make_deg(poly1, req_deg)
    poly2 = make_deg(poly2, req_deg)
    return [poly1, poly2]

def next_power_of_2(n : int):
    if n <= 0:
        return 0
    return pow(2, math.ceil(math.log(n)/math.log(2)))

def get_omega(n : int):
    assert(n > 0)
    ind = (2 * np.pi * complex(0, 1)) / n
    return np.exp(ind)
