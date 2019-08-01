
#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifdef __cplusplus
extern "C" {
#endif

#define C  299792458
#define PI 3.14159265358979323846

#define SINGLE_PRECISION 1
    
#if SINGLE_PRECISION
    #define PREC float
    #define PREC_MAX FLT_MAX
    #define PREC_SIN(x) sinf(x)
    #define PREC_COS(x) cosf(x)
    #define PREC_SQRT(x) sqrtf(x)
    #define PREC_POW(x, y) powf(x, y)
    #define PREC_ABS(x) fabsf(x)
    #define PREC_ROUND(x) roundf(x)
#else
    #define PREC double
    #define PREC_MAX DBL_MAX
    #define PREC_SIN(x) sin(x)
    #define PREC_COS(x) cos(x)
    #define PREC_SQRT(x) sqrt(x)
    #define PREC_POW(x, y) pow(x, y)
    #define PREC_ABS(x) fabs(x)
    #define PREC_ROUND(x) round(x)
#endif
    
#ifdef __cplusplus
}
#endif

#endif /* CONSTANTS_H */

