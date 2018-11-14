#ifndef __AF_BONDS_BENDING_CUH__
#define __AF_BONDS_BENDING_CUH__

#include <tuple>

#include <cuda.h>

#include "../../utilities/vector_type_helpers.cuh"

namespace af 
{
    typedef tuple<float3,float3,float3> float3_triple

    struct BendingBase
    {
        __host__ __device__
        virtual float3_triple operator()(const float3 &r, const float3 &r_behind, const float3 &r_ahead) = 0;
    };

    struct CosineBending : public BendingBase
    {
        const float K;

        __host__ __device__
        CosineBending(const float k_bend)
            : K(k_bend){}

        __host__ __device__
        float3_triple operator()(const float3 &r, const float3 &r_behind, const float3 &r_ahead)
        {
            // TODO: apply BCs to these somehow
            float3 r12 = r - r_ahead;
            float3 r23 = r_behind - r;

            float dot_r12_r23 = dot(r12, r23);
            float r12r12 = dot(r12, r12);
            float r23r23 = dot(r23, r23);
            float mag12inv = 1.0f / mag(r12);
            float mag23inv = 1.0f / mag(r23);
            float a = this->K * mag12inv * mag23inv;
            float3 A = make_float3(a, a, a);

            float3 f1 = -A * (r23 - ((dot_r12_r23 / r12r12) * r12));
            float3 f2 = -A * (((dot_r12_r23 / r12r12) * r12) - ((dot_r12_r23 / r23r23) * r23) + r12 - r23);
            float3 f3 = -A * (((dot_r12_r23 / r23r23) * r23) - r12);
            return make_tuple(f1,f2,f3);

            //return zero_float3();
        }
    };
}

#endif