#ifndef __AF_FORCES_FILAMENTKERNEL_CUH__
#define __AF_FORCES_FILAMENTKERNEL_CUH__

#include <cuda.h>
#include <device_launch_parameters.h>
#include <math_functions.h>

#include "../neighbor_finding/cells.cuh"
#include "bonds.h"
#include "pair_wise.h"
#include "../particle.cuh"

namespace af
{
    // launch a thread per filament
    __global__
    void filament_kernel(
        Particle* particles, 
        size_t n_particles,
        uint* filament_headlist,
        size_t n_filaments)
    {
        int fidx = blockIdx.x*blockDim.x + threadIdx.x;
        if (fidx >= n_filaments) return;

        // iterate over particles in filament
        uint head_particle_idx = filament_headlist[fid];
        for (Particle* part = &particles[head_particle_idx];
             part->next_idx != -1;
             part =  &particles[part->next_idx])
        {
            if (fidx == 0) printf("%d ", part->id);
        }
    }
}