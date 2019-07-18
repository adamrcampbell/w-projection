
#ifndef UTILITY_H
#define UTILITY_H

#ifdef __cplusplus
extern "C" {
#endif

    #define MIN(a,b) (((a)<(b))?(a):(b))
    #define MAX(a,b) (((a)>(b))?(a):(b))
    
    size_t get_total_ram_capacity();
    unsigned int get_next_pow_2(unsigned int x);

#ifdef __cplusplus
}
#endif

#endif /* UTILITY_H */

