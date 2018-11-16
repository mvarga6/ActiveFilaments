#ifndef __AF_FORCES_OPERATOR_H__
#define __AF_FORCES_OPERATOR_H__

#include <cuda_runtime.h>

#include "../particle.cuh"
#include "../neighbor_finding/neighbor_finder.cuh"
#include "force_kernel.cuh"
#include "preforce_kernel.cuh"

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
        void update(ParticleDeviceArray& particles, uint parts_per_filament, uint num_filaments)
        {
            const int TPB = 512;

            // update neighbors
            if (neighbors != NULL)
                neighbors->update(particles, num_filaments);

            // calculate forces
            Cells cells = neighbors->get_cells(); //TODO get this from sim container
            size_t n_cells = cells.count();

            preforce_kernel<<<TPB,n_cells/TPB + 1>>>
            (
                thrust::raw_pointer_cast(&particles[0]),
                particles.size(), 
                parts_per_filament,
                thrust::raw_pointer_cast(&cell_head_idx[0]),
                thrust::raw_pointer_cast(&cell_count[0]),
                cells, 2
            );

            cudaDeviceSynchronize();

            force_kernel<<<TPB,n_cells/TPB + 1>>>
            (
                thrust::raw_pointer_cast(&particles[0]),
                particles.size(),
                thrust::raw_pointer_cast(&cell_head_idx[0]),
                thrust::raw_pointer_cast(&cell_count[0]),
                cells, 2, opts
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