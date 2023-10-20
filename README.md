# Multiplying-polynomials

## Description
This is a simple realization of polynomial multiplication with a FFT approach.

<img src="polynomial_multiplication.pbm" alt="polynomial_multiplication.pbm" width="400"/>

The `fft.fft()` function from the `fft` module takes a poly $p ~ : ~ \deg{p} = n-1$ , a complex number $\omega = e^{\frac{2 \pi i}{n}}$ and return an array of values $p(\omega^j), ~ j \in 0 \dots n-1$.

The `fft.mult()` function multiplies two polynomials using the function described above.

## Tests
There is a `test()` function that compares the output from `fft.mult(p1, p2)` with `p1 * p2`, where `*` is the operator from `numpy.Polynomial`. `p1` and `p2` are random polynomials with integer coefficients. $1 \leqslant \deg p_1, \deg p_2 < 10$.
