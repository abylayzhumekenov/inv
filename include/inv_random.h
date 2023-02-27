#ifndef INV_RANDOM_H
#define INV_RANDOM_H

#include <pcg_variants.h>


/**
 * @brief Create a random number generator
 * 
 * @param initstate 
 * @param initseq 
 * @param unique 
 * @return pcg64_random_t 
 */
pcg64_random_t InvRandomCreate(pcg128_t initstate, pcg128_t initseq, int unique);


/**
 * @brief Generate one standard normal random number using Box-Muller transform
 * 
 * @param rng 
 * @return double 
 */
double InvRandomStdNormal(pcg64_random_t* rng);


#endif