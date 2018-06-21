#ifndef __AF_FORCES_FORCEKERNEL_CUH__
#define __AF_FORCES_FORCEKERNEL_CUH__

#include <cuda.h>
#include <device_launch_parameters.h>
#include <math_functions.h>

#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/for_each.h>


#include "../neighbor_finding/cells.cuh"
#include "bonds.h"
#include "pair_wise.h"
#include "../particle.cuh"
#include "options.cuh"

namespace af
{

    // launch a thread per cell
    __global__ void force_kernel(
        Particle* particles, 
        size_t n_particles,
        uint* cell_heads, 
        uint* cell_counts,
        Cells cells,
        uint dimension,
        ForceKernelOptions opts)
        {
            int cell_idx = blockIdx.x*blockDim.x + threadIdx.x;
            if (cell_idx < cells.count())
            {
                // create the force functors
                BondBase * backbone_bond 
                    = BondFactory::create(
                        opts.backbone_bonds,
                        opts.backbone_energy_scale,
                        opts.backbone_length_scale);

                PairWiseBase * particle_particle
                    = PairWiseFactory::create(
                        opts.filament_interaction,
                        opts.interaction_energy_scale,
                        opts.interaction_length_scale);

                BendingBase * filament_bending
                    = BondFactory::create(
                        opts.filament_bending,
                        opts.bending_energy_scale);

                // find particles for this cell
                uint head = cell_heads[cell_idx];
                uint count = cell_counts[cell_idx];

                // how many neighboring cells are there
                uint n_dirs = powf(3, dimension);

                for (uint p1_idx = head;   // the particles we're applying
                    p1_idx < head + count; // forces too
                    p1_idx++)
                {
                    // get first particle ref
                    Particle* p1 = &particles[p1_idx];
                    const float3 r1 = p1->r;

                    // aggregate forces in ths object
                    float3 f = zero_float3();

                    // positions of particles ahead and
                    // behind for calculating bending forces.
                    float3 chain_neighbor_ahead_r;
                    float3 chain_neighbor_behind_r;

                    // loop over cells
                    for (int dir = 0; dir < n_dirs; dir++)
                    {
                        // find particles for the search cell
                        uint neighbor_cell_idx = cells.neighbor_idx(cell_idx, dir);
                        if (neighbor_cell_idx == NO_CELL_IDX) continue;
                        uint neighbor_head = cell_heads[neighbor_cell_idx];
                        uint neighbor_count = cell_counts[neighbor_cell_idx];
                    
                        for (uint p2_idx = neighbor_head;        // the neighbor we're calculating
                            p2_idx < neighbor_head + neighbor_count; // forces with
                            p2_idx++)
                        {
                            // don't calculate forces with yourself
                            if (p1_idx == p2_idx) continue;

                            // get the ref to other particle
                            Particle* p2 = &particles[p2_idx];
                            const float3 r2 = p2->r;

                            // Same filament forces
                            bool same_filament = p1->filament_id == p2->filament_id;
                            if (same_filament)
                            {
                                // distance along backbone
                                int local_dist = p2->local_id - p1->local_id;
                                bool bonded = abs(local_dist) == 1;

                                if (bonded) // particles are next to eachother in chain
                                {
                                    // apply backbone bonding force
                                    f += (*backbone_bond)(r1, r2);

                                    // store positions of neighbors in chain
                                    if (local_dist == 1) // particle ahead
                                        chain_neighbor_ahead_r = r2;
                                    else 
                                        chain_neighbor_behind_r = r2;
                                } 
                                else if (local_dist > 2) // normal interaction
                                    f += (*particle_particle)(p1->r, p2->r);
                            }
                            else 
                            {
                                f += (*particle_particle)(r1, r2);

                                // TODO: extensile forces
                            }
                        }
                    }

                    // Bond Bending Forces
                    f += (*filament_bending)(r1, chain_neighbor_behind_r, chain_neighbor_ahead_r);

                    // apply forces to particle
                    p1->f += f;
                    //printf("%f %f %f\n", f.x, f.y, f.z);
                }

                delete backbone_bond;
                delete particle_particle;
                delete filament_bending;
            }
        }


    // launch a thread per cell
    __global__ void neighbor_force_kernel_by_particle(
        Particle* particles, 
        size_t n_particles,
        uint* cell_heads, 
        uint* cell_counts,
        size_t n_cells,
        uint dimension)
        {
            //int particle_idx = blockIdx.x*blockDim.x + threadIdx.x;
        }
}

#endif