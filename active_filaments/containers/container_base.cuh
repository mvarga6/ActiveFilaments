#ifndef __AF_CONTAINERS_BASE_CUH__
#define __AF_CONTAINERS_BASE_CUH__

#include <cuda.h>
#include <cuda_runtime.h>

#include "../neighbor_finding/cells.cuh"
#include "../particle.cuh"

namespace af
{
    struct ContainerBase
    {
        __host__ __device__
        ContainerBase(){}
        
        __host__ __device__
        virtual void generate_cells() = 0;

        __host__ __device__ 
        virtual void add_forces_to(Particle& particle) = 0;

        __host__ __device__
        virtual void boundary_conditions(Particle& particle) = 0;
    };
}

#endif