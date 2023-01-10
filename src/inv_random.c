#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "inv_random.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif


pcg64_random_t InvRandomCreate(pcg128_t initstate, pcg128_t initseq, int unique){

    pcg64_random_t rng;
    if(unique) pcg64_srandom_r(&rng, time(NULL), (intptr_t)&rng);  // completely unpredictable seeds
    else pcg64_srandom_r(&rng, initstate, initseq);

    return rng;
}


double InvRandomStdNormal(pcg64_random_t* rng){
    
    double u1 = ldexp(pcg64_random_r(rng), -64);
    double u2 = ldexp(pcg64_random_r(rng), -64);
    
    return sqrt(-2.0*log(u1)) * cos(2*M_PI*u2);
}
