import numpy as np
import numpy.polynomial.polynomial as polynomial
import fft

def test():
    for i in range(1, 10):
        for j in range(1, 10):
            is_eq, poly1, poly2, fft_poly, target_poly = test_multiplication(i, j)
            if not is_eq:
                print("Error!")
                print("poly1 = ", poly1)
                print("poly2 = ", poly2)
                print("fft_poly = ", fft_poly)
                print("np_poly = ", target_poly)
                return
    print("All cases passed successfully")

def test_multiplication(deg1, deg2):
    low = -10
    high = 10
    poly1 = random_poly(deg1, low, high)
    poly2 = random_poly(deg2, low, high)

    fft_poly = fft.to_round_poly(fft.mult(poly1.copy(), poly2.copy()))
    target_poly = poly1 * poly2
    [fft_poly, target_poly] = fft.align_polynomials(fft_poly, target_poly)
    is_eq = (fft_poly == target_poly)
    return [is_eq, poly1, poly2, fft_poly, target_poly]

def random_poly(deg, low, high):
    coefs = np.random.randint(low, high, deg)
    while coefs[deg - 1] == 0:
        coefs[deg - 1] = np.random.randint(low, high, 1)
    return polynomial.Polynomial(coefs)