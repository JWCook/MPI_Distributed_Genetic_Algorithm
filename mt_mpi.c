/* ========================================================================== */
/* A implementation of the Mersenne Twister algorithm for fast generation     */
/* of pseudorandom numbers, modified to support multiple parallel streams.    */
/* It returns random random integers in the range 0 to 2^32-1 with a period   */
/* of2^19937-1.                                                               */
/* ========================================================================== */

/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "mt_mpi.h"


static unsigned long mt[NUM_RNG][N]; // State vector
static int mti[NUM_RNG];       // mti[x]==N+1 means mt[x][N] is not initialized


/* -------------------------------------------------------------------------- */
/* Generate a seed value                                                      */
/* -------------------------------------------------------------------------- */
unsigned long time_seed() {
    time_t cur_time = time (0);
    unsigned char *p = (unsigned char *)&cur_time;
    unsigned long seed = 0;
    size_t i;

    for (i = 0; i < sizeof cur_time; i++)
        seed = seed * (UCHAR_MAX + 2U) + p[i];

    return seed;
}


/* -------------------------------------------------------------------------- */
/* Initialize state vector                                                    */
/* my_rank        : The processor to initialize a for                         */
/* seed            : A seed value to use; generate a new seed if this is 0.   */
/* -------------------------------------------------------------------------- */
void mt_init(int my_rank, unsigned long seed) {
    int i;

    if (seed == 0) seed = time_seed();
    seed += my_rank;
    for (i = 0; i < NUM_RNG; i++) mti[i] = N+1;
    mt[my_rank][0]= seed & 0xffffffffU;
    for (mti[my_rank] = 1; mti[my_rank] < N; mti[my_rank]++)
        mt[my_rank][mti[my_rank]] = (69069 * mt[my_rank][mti[my_rank]-1]) &
        0xffffffffU;
}


/* -------------------------------------------------------------------------- */
/* Generate a random integer on the interval [0,0xffffffff]                   */
/* my_rank        : The processor to generate a random number for             */
/* -------------------------------------------------------------------------- */
unsigned long mt_rand_32(int my_rank) {
    int    i, kk;
    unsigned long y;

    static int first = 1;
    static unsigned long mag01[NUM_RNG][2];

    // Init mag01
    if (first == 1) {
        for (i = 0; i < NUM_RNG; i++) {
            mag01[i][0] = 0x0;
            mag01[i][1] = MATRIX_A;
        }
        first = 0;
    }

    // Generate N words at a time
    if (mti[my_rank] >= N) {
        if (mti[my_rank] == N+1) mt_init(my_rank, 0);

        for (kk=0;kk<N-M;kk++) {
            y = (mt[my_rank][kk]&UPPER_MASK)|(mt[my_rank][kk+1]&LOWER_MASK);
            mt[my_rank][kk] = mt[my_rank][kk+M] ^
                (y >> 1) ^ mag01[my_rank][y & 0x1];
        }

        for (;kk<N-1;kk++) {
            y = (mt[my_rank][kk]&UPPER_MASK)|(mt[my_rank][kk+1]&LOWER_MASK);
            mt[my_rank][kk] = mt[my_rank][kk+(M-N)] ^
                (y >> 1) ^ mag01[my_rank][y & 0x1];
        }

        y = (mt[my_rank][N-1]&UPPER_MASK)|(mt[my_rank][0]&LOWER_MASK);
        mt[my_rank][N-1] = mt[my_rank][M-1] ^
                (y >> 1) ^ mag01[my_rank][y & 0x1];
        mti[my_rank] = 0;
    }

    y = mt[my_rank][mti[my_rank]++];
    y ^= (y >> 11);
    y ^= (y << 7)  & TEMPERING_MASK_B;
    y ^= (y << 15) & TEMPERING_MASK_C;
    y ^= (y >> 18);

    return y;
}


/* -------------------------------------------------------------------------- */
/* Generate a random real number on the interval [0,n]                        */
/* -------------------------------------------------------------------------- */
long double mt_rand_real(double n, int my_rank) {
    return mt_rand_32(my_rank)*(n/MT_MAX);
}


/* -------------------------------------------------------------------------- */
/* Generate a random integer on the interval [0,n]                            */
/* -------------------------------------------------------------------------- */
int mt_rand_int(int n, int my_rank) {
    double full_value = mt_rand_32(my_rank)*(n/MT_MAX);
    return rint(full_value);
}


/* -------------------------------------------------------------------------- */
/* Generate a random bit                                                      */
/* -------------------------------------------------------------------------- */
char mt_rand_bit(int my_rank) {
    double full_value = mt_rand_32(my_rank)*(1.0/MT_MAX);
    int rounded = rint(full_value);
    return rounded + '0';
}


/* -------------------------------------------------------------------------- */
/* Return true with the specified probability                                 */
/* -------------------------------------------------------------------------- */
int mt_probability(float prob, int my_rank) {
    return (mt_rand_real(1, my_rank) <= prob);
}

