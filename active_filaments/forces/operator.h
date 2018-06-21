#ifndef __AF_FORCES_OPERATOR_H__
#define __AF_FORCES_OPERATOR_H__

#include <cuda_runtime.h>

#include "../particle.cuh"
#include "../neighbor_finding/neighbor_finder.cuh"
#include "force_kernel.cuh"

namespace af
{
    class Forces
    {
        NeighborFinder* neighbors;
        ForceKernelOptions opts;
    public:
        __host__
        Forces(ForceKernelOptions options, NeighborFinder* neighbor_finder = NULL)
            : neighbors(neighbor_finder), opts(options){}

        __host__ 
        void update(ParticleDeviceArray& particles)
        {
            // update neighbors
            if (neighbors != NULL)
                neighbors->update(particles);

            // calculate forces
            Cells cells = neighbors->get_cells();
            size_t n_cells = cells.count();

            force_kernel<<<256,n_cells/256 + 1>>>
            (
                thrust::raw_pointer_cast(&particles[0]),
                particles.size(),
                thrust::raw_pointer_cast(&cell_head_idx[0]),
                thrust::raw_pointer_cast(&cell_count[0]),
                cells, 3, opts
            );
        }

        __host__ 
        void set_options(ForceKernelOptions new_options)
        {
            opts = new_options;
        }
    };
}

#endif