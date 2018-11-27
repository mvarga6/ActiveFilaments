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
    //
    // Returns pointer to next particle in filament.
    // Return NULL if particle is tail.
    //
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
        size_t n_filaments,
        ForceKernelOptions opts)
    {
        int fidx = blockIdx.x*blockDim.x + threadIdx.x;
        if (fidx >= n_filaments) return;

        // create the force functors
        BondBase * backbone_bond 
            = BondFactory::create(
                opts.backbone_bonds,
                opts.backbone_energy_scale,
                opts.backbone_length_scale);

        BendingBase * filament_bending
            = BondFactory::create(
                opts.filament_bending,
                opts.bending_energy_scale);

        // Gather particle pointers for this filament
        thrust::device_vector<Particle*> filament;
        uint head_particle_idx = filament_headlist[fidx];
        for (Particle* p = &particles[head_particle_idx];
            p != NULL; p = next(particles, p))
            filament.push_back(p);

        // filament lenght (support filament of difference sizes)
        uint n = filament.size();

        // Bending forces
        for (int i = 0; i+2 < n; i++)
            (*filament_bending)(filament[i], filament[i+1], filament[i+2]);

        // Bonding forces
        for (int i = 0; i+1 < n; i++)
            (*backbone_bond)(filament[i], filament[i+1]);

        //
        // TODO: DRIVE ACTIVITY
        //

        // Cleanup local memory
        delete backbone_bond;
        delete filament_bending;
    }
}

#endif