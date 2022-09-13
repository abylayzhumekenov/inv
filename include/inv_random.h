#ifndef INV_RANDOM_H
#define INV_RANDOM_H

#include <pcg_variants.h>


/**
 * @brief Generate one standard normal random number using Box-Muller transform
 * 
 * @param rng 
 * @return double 
 */
double random_std_normal(pcg64_random_t* rng);


#endif