#ifndef __AF_FORCES_PAIRWISE_BASE_CUH__
#define __AF_FORCES_PAIRWISE_BASE_CUH__

#include <cuda.h>

namespace af 
{
    struct PairWiseBase
    {
        __host__ __device__
        virtual float3 operator()(const float3& r1, const float3& r2) = 0;
    };

    struct NoInteraction : public PairWiseBase
    {
        __host__ __device__
        NoInteraction(){}

        __host__ __device__
        float3 operator()(const float3 &r1, const float3 &r2)
        {
            return zero_float3(); 
        }
    };
}

#endif