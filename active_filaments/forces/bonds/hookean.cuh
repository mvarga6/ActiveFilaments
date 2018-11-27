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
        float3 operator()(Particle* p1, Particle* p2)
        {
            // TODO: Boundary Conditions
            float3 dr = p2->r - p1->r;
            const float r = mag(dr);
            const float f = -K * (R0 - r) / r;
            p1.f += f * dr;
            p2.f -= f * dr;
        }
    };
}

#endif