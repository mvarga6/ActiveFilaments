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
    __device__
    Particle* next(Particle* particles, Particle* current)
    {
        if (current->next_idx == -1)
            return NULL;
        else
            return &particles[current->next_idx];
    }

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
        if (fidx == 0)
        {
            uint head_particle_idx = filament_headlist[fidx];
            for (Particle* p = &particles[head_particle_idx];
                p != NULL; p = next(particles, p))
            {
                printf("%d ", p->id);
            }
            printf("\n");
        }
    }
}

#endif