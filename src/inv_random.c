#include <math.h>
#include "inv_random.h"


double random_std_normal(pcg64_random_t* rng){
    
    double u1 = ldexp(pcg64_random_r(rng), -64);
    double u2 = ldexp(pcg64_random_r(rng), -64);
    
    return sqrt(-2.0*log(u1)) * cos(2*M_PI*u2);
}