
#include <math.h>
#include <float.h>
#include <stdbool.h>

#include "complex.h"
#include "constants.h"

Complex complex_add(Complex z1, Complex z2)
{
    return (Complex) {
        .real = z1.real + z2.real,
        .imag = z1.imag + z2.imag
    };
}

Complex complex_subtract(Complex z1, Complex z2)
{
    return (Complex) {
        .real = z1.real - z2.real,
        .imag = z1.imag - z2.imag
    };
}

Complex complex_multiply(Complex z1, Complex z2)
{
    return (Complex) {
        .real = z1.real * z2.real - z1.imag * z2.imag,
        .imag = z1.imag * z2.real + z1.real * z2.imag
    };
}

Complex complex_scale(Complex z, PREC scalar)
{
    return (Complex) {
        .real = z.real * scalar,
        .imag = z.imag * scalar
    };
}

Complex complex_conj_exp(PREC phase)
{
    return (Complex) {
        .real = PREC_COS(2.0 * PI * phase),
        .imag = -PREC_SIN(2.0 * PI * phase)
    };
}

PREC complex_magnitude(Complex z)
{
    return PREC_SQRT(z.real * z.real + z.imag * z.imag);
}