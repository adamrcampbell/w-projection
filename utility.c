
#include <stdlib.h>
#include <stdbool.h>
#include <sys/sysinfo.h> // Linux specific

size_t get_total_ram_capacity()
{
    struct sysinfo memInfo;
    sysinfo(&memInfo);
    return memInfo.totalram * memInfo.mem_unit;
}

// Suitable for 32-bit unsigned integers
unsigned int get_next_pow_2(unsigned int x)
{
    x--;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x++;
    
    return x;
}

bool is_power_of_two(unsigned int x)
{
    return ((x & ~(x-1))==x)? x : 0;
}
