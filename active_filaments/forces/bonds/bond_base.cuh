#ifndef __AF_FORCES_BOND_BASE_CUH__
#define __AF_FORCES_BOND_BASE_CUH__

#include <cuda.h>

#include "../particle.cuh"

namespace af
{
    struct BondBase
    {
        __host__ __device__
        virtual void operator()(Particle* p1, Particle* p2) = 0;
    };

    struct NoBond : public BondBase
    {
        __host__ __device__
        NoBond(){}

        __host__ __device__
        void operator()(Particle* p1, Particle* p2)
        {
            
        }
    };
}

#endif