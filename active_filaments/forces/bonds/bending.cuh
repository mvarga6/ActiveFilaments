#ifndef __AF_BONDS_BENDING_CUH__
#define __AF_BONDS_BENDING_CUH__


#include <cuda.h>

#include "../../utilities/vector_type_helpers.cuh"
#include "../particle.cuh"

namespace af 
{
    struct BendingBase
    {
        __host__ __device__
        virtual void operator()(Particle *p1, Particle *p2, Particle *p3) = 0;
    };

    struct CosineBending : public BendingBase
    {
        const float K;

        __host__ __device__
        CosineBending(const float k_bend)
            : K(k_bend){}

        __host__ __device__
        void operator()(Particle *p1, Particle *p2, Particle *p3)
        {
            // TODO: apply BCs to these somehow
            float3 r12 = p2->r - p1->r;
            float3 r23 = p3->r - p2->r;

            float dot_r12_r23 = dot(r12, r23);
            float r12r12 = dot(r12, r12);
            float r23r23 = dot(r23, r23);
            float mag12inv = 1.0f / mag(r12);
            float mag23inv = 1.0f / mag(r23);
            float A = -this->K * mag12inv * mag23inv;

            float c1 = dot_r12_r23 / r12r12;
            float c2 = dot_r12_r23 / r23r23;

            p1->f += A * (r23 - c1*r12);
            p2->f += A * ((c1*r12 + r12) - (c2*r23 + r23));
            p3->f += A * (c2*r23 - r12);
        }
    };
}

#endif