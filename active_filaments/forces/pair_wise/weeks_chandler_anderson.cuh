#ifndef __AF_FORCES_WEEKSCHANDLERANDERSON_CUH__
#define __AF_FORCES_WEEKSCHANDLERANDERSON_CUH__

#include <cuda.h>

#include "../../particle.cuh"
#include "lennard_jones.cuh"

namespace af 
{
    struct WeeksChandlerAnderson : public LennardJones
    {
        __host__ __device__
        WeeksChandlerAnderson(float sigma, float epsilon)
            : LennardJones(sigma, epsilon)
            {
                // change the cutoff
                const float rmin = powf(2.0, 1.0f/6.0f) * sigma;
                _rrcut = rmin * rmin;
            }

            // Don't need to override operator() becauce
            // we've already set the correct cutoff.
            // All the behavior is the same besides that.
    };
}

#endif