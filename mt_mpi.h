#ifndef MT_LF_H_
#define MT_LF_H

unsigned long    time_seed(void);
void             mt_init(int, unsigned long);
unsigned long    mt_rand_32(int);
long double      mt_rand_real(double, int);
int              mt_rand_int(int, int);
char             mt_rand_bit(int);
int              mt_probability(float, int);

#endif


#define NUM_RNG             128
#define N                   624
#define M                   397
#define MT_MAX              4294967295.0
#define TEMPERING_MASK_B    0x9d2c5680U
#define TEMPERING_MASK_C    0xefc60000U
#define MATRIX_A            0x9908b0dfU        // Constant vector a
#define UPPER_MASK          0x80000000U        // Most significant w-r bits
#define LOWER_MASK          0x7fffffffU        // Least significant r bits
