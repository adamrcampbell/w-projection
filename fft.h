
#ifndef FFT_H
#define FFT_H

#ifdef __cplusplus
extern "C" {
#endif

    void fft_2d(Complex *matrix, int number_channels);
    void fft_shift_2d(Complex *matrix, int size);
    void calc_bit_reverse_indices(int n, int* indices);
    
    void inverseFFT2dVectorRadixTransform(int numChannels, Complex *input, Complex *output);
    void calcBitReversedIndices(int n, int* indices);
    void fft2dShift(int n, Complex *input, Complex *shifted);

#ifdef __cplusplus
}
#endif

#endif /* FFT_H */

