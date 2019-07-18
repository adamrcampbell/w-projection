
#ifndef W_PROJECTION_H
#define W_PROJECTION_H

#ifdef __cplusplus
extern "C" {
#endif
   
    void generate_w_projection_kernels(void);
    void generate_phase_screen(int iw, int conv_size, int inner, PREC sampling, PREC w_scale, PREC* taper, Complex *screen);
    void normalize_kernels_by_maximum(Complex *kernels, PREC *maximums, int number_w_planes, int conv_half_size);
    void normalize_kernels_sum_of_one(Complex *kernels, int number_w_planes, int conv_half_size, int oversample);
    
    void normalize_plane(Complex *plane, int resolution);
    void crop_plane(Complex *plane, Complex *cropped_plane, int resolution, int support);
    PREC calculate_support(PREC w, int min_support, PREC w_max_support_ratio);
    void save_kernel_to_file(Complex *plane, int resolution, int clipping, char *real, char *imag);
    
#ifdef __cplusplus
}
#endif

#endif /* W_PROJECTION_H */

