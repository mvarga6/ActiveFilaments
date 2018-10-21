#ifndef __AF_BONDS_HOOKEAN_CUH__
#define __AF_BONDS_HOOKEAN_CUH__

#include <cuda.h>

#include "bond_base.cuh"
#include "../../utilities/vector_type_helpers.cuh"

namespace af
{
    struct HookeanBond : public BondBase
    {
        const float K, R0;

        __host__ __device__
        HookeanBond(float k, float r0)
            : K(k), R0(r0){}

        __host__ __device__
        float3 operator()(const float3 &r1, const float3 &r2)
        {
            float3 dr = r2 - r1;
            const float r = mag(dr);
            const float f = -K * (R0 - r) / r;
            return dr * f;
        }
    };
}

#endif