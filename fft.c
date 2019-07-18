
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "constants.h"
#include "complex.h"
#include "fft.h"

void fft_2d(Complex *matrix, int number_channels)
{
    // Calculate bit reversed indices
    int* bit_reverse_indices = calloc(number_channels, sizeof(int));
    calc_bit_reverse_indices(number_channels, bit_reverse_indices);
    Complex *reverse_buffer = calloc(number_channels * number_channels, sizeof(Complex));
    
    for(int row = 0; row < number_channels; ++row)
        for(int col = 0; col < number_channels; ++col)
        {
            int row_reverse = bit_reverse_indices[row];
            int col_reverse = bit_reverse_indices[col];
            int bit_reverse_index = row_reverse * number_channels + col_reverse;
            int matrix_index = row * number_channels + col;
            reverse_buffer[matrix_index] = matrix[bit_reverse_index];
            // printf("%d -> %d\n", matrix_index, bit_reverse_index);
        }
    
    memcpy(matrix, reverse_buffer, number_channels * number_channels * sizeof(Complex));
    free(reverse_buffer);
    free(bit_reverse_indices);
    
    for(int m = 2; m <= number_channels; m *= 2)
    {
        Complex omegaM = (Complex) {.real = PREC_COS(PI * 2.0 / m), .imag = PREC_SIN(PI * 2.0 / m)};
        
        for(int k = 0; k < number_channels; k += m)
        {
            for(int l = 0; l < number_channels; l += m)
            {
                Complex x = (Complex) {.real = 1.0, .imag = 0.0};
                
                for(int i = 0; i < m / 2; i++)
                {
                    Complex y = (Complex) {.real = 1.0, .imag = 0.0};
                    
                    for(int j = 0; j < m / 2; j++)
                    {   
                        // Perform 2D butterfly operation in-place at (k+j, l+j)
                        int in00Index = (k+i) * number_channels + (l+j);
                        Complex in00 = matrix[in00Index];
                        int in01Index = (k+i) * number_channels + (l+j+m/2);
                        Complex in01 = complex_multiply(matrix[in01Index], y);
                        int in10Index = (k+i+m/2) * number_channels + (l+j);
                        Complex in10 = complex_multiply(matrix[in10Index], x);
                        int in11Index = (k+i+m/2) * number_channels + (l+j+m/2);
                        Complex in11 = complex_multiply(complex_multiply(matrix[in11Index], x), y);
                        
                        Complex temp00 = complex_add(in00, in01);
                        Complex temp01 = complex_subtract(in00, in01);
                        Complex temp10 = complex_add(in10, in11);
                        Complex temp11 = complex_subtract(in10, in11);
                        
                        matrix[in00Index] = complex_add(temp00, temp10);
                        matrix[in01Index] = complex_add(temp01, temp11);
                        matrix[in10Index] = complex_subtract(temp00, temp10);
                        matrix[in11Index] = complex_subtract(temp01, temp11);
                        y = complex_multiply(y, omegaM);
                    }
                    x = complex_multiply(x, omegaM);
                }
            }
        }
    }
    
    for(int row = 0; row < number_channels; ++row)
        for(int col = 0; col < number_channels; ++col)
        {   
            int matrix_index = row * number_channels + col;
            PREC reciprocal = 1.0 / (number_channels * number_channels);
            matrix[matrix_index] = complex_scale(matrix[matrix_index], reciprocal);
        }
}

void calc_bit_reverse_indices(int n, int* indices)
{   
    for(int i = 0; i < n; ++i)
    {
        // Calculate index r to which i will be moved
        unsigned int i_prime = i;
        int r = 0;
        for(int j = 1; j < n; j *= 2)
        {
            int b = i_prime & 1;
            r = (r << 1) + b;
            i_prime = (i_prime >> 1);
        }
        indices[i] = r;
    }
}

void fft_shift_2d(Complex *matrix, int size)
{
    for(int row = 0; row < size; ++row)
        for(int col = 0; col < size; ++col)
        {
            int matrix_index = row * size + col;
            PREC scalar = 1 - 2 * ((row + col) & 1);
            matrix[matrix_index] = complex_scale(matrix[matrix_index], scalar);
        }
}
