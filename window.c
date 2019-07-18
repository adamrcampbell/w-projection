
#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "window.h"

void populate_ps_window(PREC *window, int size)
{
    for(int index = 0; index < size; ++index)
    {
        PREC nu = PREC_ABS(calculate_window_stride(index, size));
        window[index] = prolate_spheroidal(nu);
    }
}

PREC calculate_window_stride(int index, int size)
{
    return (index - size / 2) / ((PREC) size / 2.0);
}

PREC prolate_spheroidal(PREC nu)
{   
    static PREC p[] = {0.08203343, -0.3644705, 0.627866, -0.5335581, 0.2312756,
        0.004028559, -0.03697768, 0.1021332, -0.1201436, 0.06412774};
    static PREC q[] = {1.0, 0.8212018, 0.2078043,
        1.0, 0.9599102, 0.2918724};
    
    int part, sp, sq;
    PREC nuend, delta, top, bottom;
    
    if(nu >= 0.0 && nu < 0.75)
    {
        part = 0;   
        nuend = 0.75;
    }
    else if(nu >= 0.75 && nu <= 1.0)
    {
        part = 1;
        nuend = 1.0;
    }
    else
        return 0.0;

    delta = nu * nu - nuend * nuend;
    sp = part * 5;
    sq = part * 3;
    top = p[sp];
    bottom = q[sq];
    
    for(int i = 1; i < 5; i++)
        top += p[sp+i] * PREC_POW(delta, i);
    for(int i = 1; i < 3; i++)
        bottom += q[sq+i] * PREC_POW(delta, i);
    return (bottom == 0.0) ? 0.0 : top/bottom;
}
