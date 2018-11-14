#ifndef __AF_FORCES_BOND_BASE_CUH__
#define __AF_FORCES_BOND_BASE_CUH__

#include <cuda.h>

namespace af
{
    struct BondBase
    {
        __host__ __device__
        virtual float3 operator()(const float3& r1, const float3& r2) = 0;
    };

    struct NoBond : public BondBase
    {
        __host__ __device__
        NoBond(){}

        __host__ __device__
        float3 operator()(const float3 &r1, const float3 &r2)
        {
            return zero_float3(); 
        }
    };
}

#endif