
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "complex.h"
#include "fft.h"
#include "w_projection.h"

int main(int argc, char** argv) {
    
    generate_w_projection_kernels();
    
    return (EXIT_SUCCESS);
}
