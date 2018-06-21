#ifndef __AF_FORCES_LENNARDJONES_CUH__
#define __AF_FORCES_LENNARDJONES_CUH__

#include <cuda.h>

#include "../../particle.cuh"
#include "pairwise_base.cuh"

namespace af 
{
    struct LennardJones : public PairWiseBase
    {
        float _2sigma6;
        float _A;
        float _rrcut;

        __host__ __device__
        LennardJones(float sigma, float epsilon)
        {
            _2sigma6 = 2.0f * powf(sigma, 6);
            _A = 24.0f * epsilon * _2sigma6;
            _rrcut = 2.5f * sigma;
        }

        __host__ __device__
        virtual float3 operator()(const float3& r1, const float3& r2)
        {
             float3 dr = r2 - r1;
            const float rr = mag_sqrd(dr);
            
            if (rr > _rrcut)
                return zero_float3();

            const float r8 = rr*rr*rr*rr;
            const float r14 = r8*rr*rr*rr;
            const float f = _A * ((_2sigma6 / r14) - (1.0f / r8));
            return dr * f;
        }
    };
}

#endif