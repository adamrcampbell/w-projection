
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

#include "constants.h"
#include "complex.h"
#include "fft.h"
#include "window.h"
#include "utility.h"
#include "w_projection.h"

void generate_w_projection_kernels(void)
{
    if(SINGLE_PRECISION)
        printf(">>> INFO: Generating W-Projection kernels using single precision...\n");
    else
        printf(">>> INFO: Generating W-Projection kernels using double precision...\n");
    
    int number_w_planes = 339;
    int grid_size = 18000;
    int image_size = 15000;
    
    int oversample = 4;
    int min_support = 4;
    int max_support = 44;
    
    size_t max_bytes_per_plane = 12 * 1024 * 1024; // 12MB
    
    PREC max_uvw  = 7083.386050;
    PREC w_scale = PREC_POW(number_w_planes - 1, 2.0) / max_uvw;
    PREC cell_size = 6.39954059065e-06;
    PREC w_to_max_support_ratio = (max_support - min_support) / max_uvw;
    PREC fov = cell_size * image_size;
    
    // Calculate convolution kernel memory requirements
    size_t max_mem_bytes = MIN(max_bytes_per_plane * number_w_planes, get_total_ram_capacity());
    PREC max_conv_size = PREC_SQRT(max_mem_bytes / (16.0 * number_w_planes));
    printf("Max conv size: %f\n", max_conv_size);
    int nearest = get_next_pow_2((unsigned int) 2 * (int) (max_conv_size / 2.0));
    printf("Nearest: %d\n", nearest);
    int conv_size = nearest;
    printf("Conv size: %d\n", conv_size);
    int conv_half_size = conv_size / 2;
    printf("Conv half size: %d\n", conv_half_size);
    
    int inner = conv_size / oversample;
    PREC max_l = PREC_SIN(0.5 * fov);
    PREC sampling = ((2.0 * max_l * oversample) / image_size) * ((PREC) grid_size / (PREC) conv_size);
    
    // Allocation of memory
    Complex *kernels = calloc(number_w_planes * conv_half_size * conv_half_size, sizeof(Complex));
    Complex *screen = calloc(conv_size * conv_size, sizeof(Complex));
    PREC *maximums = calloc(number_w_planes, sizeof(PREC));
    
    PREC* taper = calloc(inner, sizeof(PREC));
    populate_ps_window(taper, inner);
    
    printf(">>> UPDATE: Creating w projection kernels...\n");
    for(int iw = 0; iw < number_w_planes; ++iw)
    {
        printf(">>> UPDATE: Creating kernel %d\n", iw);
        // Zero out screen
        memset(screen, 0, conv_size * conv_size * sizeof(Complex));
        
        // Generate screen
        generate_phase_screen(iw, conv_size, inner, sampling, w_scale, taper, screen);
        
        printf(">>> UPDATE: Executing Fourier Transform...\n");
        // FFT
        fft_2d(screen, conv_size);
        
        // store maximum
        maximums[iw] = complex_magnitude(screen[0]);
        
        printf(">>> UPDATE: Clipping useful quadrant for further processing...\n");
        // Clip
        for(int row = 0; row < conv_half_size; ++row)
            for(int col = 0; col < conv_half_size; ++col)
            {
                int offset = iw * conv_half_size * conv_half_size;
                int k_index = offset + row * conv_half_size + col;
                kernels[k_index] = screen[row * conv_size + col];
            }
    }
    
    free(taper);
    free(screen);
    
    printf(">>> UPDATE: Normalizing quadrants by the global maximum (typically peak of plane where w == 0)...\n");
    normalize_kernels_by_maximum(kernels, maximums, number_w_planes, conv_half_size);
    free(maximums);
    
    printf(">>> UPDATE: Normalizing quadrants to the sum of one...\n");
    normalize_kernels_sum_of_one(kernels, number_w_planes, conv_half_size, oversample);
    
    printf(">>> UPDATE: Extracting oversampled support kernel from quadrants...\n");
    FILE *kernel_real_file = fopen("kernels/w-proj_kernels_real.csv", "w");
    FILE *kernel_imag_file = fopen("kernels/w-proj_kernels_imag.csv", "w");
    FILE *support_file = fopen("kernels/w-proj_supports.csv", "w");
    
    for(int iw = 0; iw < number_w_planes; ++iw)
    {
        PREC w = iw * iw / w_scale;
        PREC support = calculate_support(w, min_support, w_to_max_support_ratio);
        int oversampled_support = (PREC_ROUND(support) + 1) * oversample;
        int kernel_offset = iw * conv_half_size * conv_half_size;
        
        fprintf(support_file, "%d\n", (int) PREC_ROUND(support));
        
        for(int row = 0; row < oversampled_support; ++row)
        {
            for(int col = 0; col < oversampled_support; ++col)
            {
                int plane_index = kernel_offset + row * conv_half_size + col;
                
                fprintf(kernel_real_file, "%.10f ", kernels[plane_index].real);
                fprintf(kernel_imag_file, "%.10f ", kernels[plane_index].imag);
            }
        }
        
        fprintf(kernel_real_file, "\n");
        fprintf(kernel_imag_file, "\n");
    }
    
    fclose(support_file);
    fclose(kernel_imag_file);
    fclose(kernel_real_file);
    free(kernels);
    
    printf(">>> UPDATE: W-Projection kernels successfully created, exiting...\n");
}

void normalize_kernels_sum_of_one(Complex *kernels, int number_w_planes, int conv_half_size, int oversample)
{
    PREC sum = 0.0;
    
    for (int iy = -4; iy <= 4; ++iy)
        for (int ix = -4; ix <= 4; ++ix)
            sum += kernels[(abs(ix) * oversample + conv_half_size * (abs(iy) * oversample))].real;
    
    unsigned int number_of_samples = number_w_planes * conv_half_size * conv_half_size;
    for(unsigned int index = 0; index < number_of_samples; ++index)
        kernels[index] = complex_scale(kernels[index], 1.0 / sum);
}

void normalize_kernels_by_maximum(Complex *kernels, PREC *maximums, int number_w_planes, int conv_half_size)
{
    PREC maximum = -PREC_MAX;
    for (int iw = 0; iw < number_w_planes; ++iw)
    {
        maximum = MAX(maximum, maximums[iw]);
        printf(">>> Plane index: %d = %f\n", iw, maximums[iw]);
    }
    
    unsigned int number_of_samples = number_w_planes * conv_half_size * conv_half_size;
    for(unsigned int index = 0; index < number_of_samples; ++index)
        kernels[index] = complex_scale(kernels[index], 1.0 / maximum);
}

void generate_phase_screen(int iw, int conv_size, int inner, PREC sampling, PREC w_scale, PREC* taper, Complex *screen)
{
    PREC f = (2.0 * PI * iw * iw) / w_scale;
    int inner_half = inner / 2;
    
    for(int iy = -inner_half; iy < inner_half; ++iy)
    {
        PREC taper_y = taper[iy + inner_half];
        PREC m = sampling * (PREC) iy;
        PREC msq = m*m;
        int offset = (iy > -1 ? iy : (iy + conv_size)) * conv_size;
        
        for(int ix = -inner_half; ix < inner_half; ++ix)
        {
            PREC l = sampling * (PREC) ix;
            PREC rsq = l * l + msq;
            if (rsq < 1.0) {
                PREC taper_x = taper[ix + inner_half];
                PREC taper = taper_x * taper_y;
                int index = (offset + (ix > -1 ? ix : (ix + conv_size)));
                PREC phase = f * (PREC_SQRT(1.0 - rsq) - 1.0);
                screen[index] = (Complex) {
                    .real = taper * PREC_COS(phase),
                    .imag = taper * PREC_SIN(phase)
                };
            }
        }
    }
}

void crop_plane(Complex *plane, Complex *cropped_plane, int resolution, int support)
{
    int half_resolution = resolution/2;
    int half_support = support/2;
    
    for(int row_index = -half_support; row_index <= half_support; ++row_index)
    {
        int cropped_row = (row_index + half_support) * support;
        int plane_row = (row_index + half_resolution) * resolution;
        
        for(int col_index = -half_support; col_index <= half_support; ++col_index)
        {
            int cropped_index = cropped_row + half_support + col_index;
            int plane_index = plane_row + half_resolution + col_index;
            
            cropped_plane[cropped_index] = plane[plane_index];
        }
    }
}

PREC calculate_support(PREC w, int min_support, PREC w_max_support_ratio)
{
    // int support = (int) (PREC_ABS(w_max_support_ratio * w) + min_support);
    // return (support % 2 == 0) ? support + 1 : support;
    
    return PREC_ABS(w_max_support_ratio * w) + min_support;
}
