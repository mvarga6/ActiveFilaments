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
        void operator()(Particle* p1, Particle* p2)
        {
            // TODO: Boundary Conditions
            float3 dr = p2->r - p1->r;
            const float rr = mag_sqrd(dr);
            const float f = H / (1 - rr/RRmax);
            p1->f += f * dr;
            p2->f -= f * dr;
        }
    };
}

#endif