
#ifndef COMPLEX_H
#define COMPLEX_H

#ifdef __cplusplus
extern "C" {
#endif
    
    #include <stdbool.h>
    #include "constants.h"
    
    typedef struct Complex {
        PREC real;
        PREC imag;
    } Complex;
    
    Complex complex_add(Complex z1, Complex z2);
    Complex complex_subtract(Complex z1, Complex z2);
    Complex complex_multiply(Complex z1, Complex z2);
    Complex complex_scale(Complex z, PREC scalar);
    Complex complex_conj_exp(PREC phase);
    PREC complex_magnitude(Complex z);
    bool complex_approx_equal(Complex z1, Complex z2);

#ifdef __cplusplus
}
#endif

#endif /* COMPLEX_H */

