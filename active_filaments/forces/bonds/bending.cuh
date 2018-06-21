#ifndef __AF_BONDS_BENDING_CUH__
#define __AF_BONDS_BENDING_CUH__

#include <cuda.h>

namespace af 
{
    struct BendingBase
    {
        __host__ __device__
        virtual float3 operator()(const float3 &r, const float3 &r_behind, const float3 &r_ahead) = 0;
    };

    struct CosineBending : public BendingBase
    {
        const float K;

        __host__ __device__
        CosineBending(const float k_bend)
            : K(k_bend){}

        __host__ __device__
        float3 operator()(const float3 &r, const float3 &r_behind, const float3 &r_ahead)
        {
            return zero_float3();
        }
    };
}

#endif