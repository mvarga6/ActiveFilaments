#ifndef __AF_BONDS_FENE_CUH__
#define __AF_BONDS_FENE_CUH__

#include <cuda.h>

#include "bond_base.cuh"
#include "../../utilities/vector_type_helpers.cuh"

namespace af
{
    struct FeneBond : public BondBase
    {
        const float H; 
        float RRmax;

        __host__ __device__
        FeneBond(float h, float r_max)
            : H(h)
            {
                RRmax = r_max * r_max;
            }

        __host__ __device__
        float3 operator()(const float3 &r1, const float3 &r2)
        {
            float3 dr = r2 - r1;
            const float rr = mag_sqrd(dr);
            const float f = H / (1 - rr/RRmax); 
            return f * dr;
        }
    };
}

#endif