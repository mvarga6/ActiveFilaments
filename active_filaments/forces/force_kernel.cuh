#ifndef __AF_FORCES_FORCEKERNEL_CUH__
#define __AF_FORCES_FORCEKERNEL_CUH__

#include <cuda.h>
#include <device_launch_parameters.h>
#include <math_functions.h>

// #include <thrust/device_vector.h>
// #include <thrust/sequence.h>
// #include <thrust/for_each.h>

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

                // local variables so were not accessing
                // global memory for every calculation
                float3 r1, f;
                uint p1_local_id;
                int next_idx, prev_idx;
                float3 next_r, prev_r;

                for (uint p1_idx = head;   // the particles we're applying
                    p1_idx < head + count; // forces too
                    p1_idx++)
                {
                    // get first particle ref
                    Particle *p1 = &particles[p1_idx];
                    r1 = p1->r;
                    next_idx = p1->next_idx;
                    prev_idx = p1->prev_idx;
                    p1_local_id = p1->local_id;

                    // aggregate forces in ths object
                    f = zero_float3();

                    //
                    // Filament Backbone Forces
                    //
                    
                    if (next_idx >= 0) // with particle ahead
                    {
                        next_r = particles[next_idx].r;
                        f += (*backbone_bond)(r1, next_r);
                    }
                        
                    if (prev_idx >= 0) // with particle behind
                    {
                        prev_r = particles[prev_idx].r;
                        f += (*backbone_bond)(r1, prev_r);
                    }
                        
                    //
                    // Bond Bending Forces
                    // 

                    float3_triple f1f2f3 = (*filament_bending)(r1, prev_r, next_r);
                    particles[prev_idx].f += get<0>(f1f2f3);
                    f += get<1>(f1f2f3);
                    particles[next_idx].f += get<2>(f1f2f3);
  
                    //
                    // Particle-Particle Forces
                    //

                    for (int dir = 0; dir < n_dirs; dir++) // loop over cells
                    {
                        // find particles for the search cell
                        uint neighbor_cell_idx = cells.neighbor_idx(cell_idx, dir);
                        if (neighbor_cell_idx == NO_CELL_IDX) continue;
                        uint neighbor_head = cell_heads[neighbor_cell_idx];
                        uint neighbor_count = cell_counts[neighbor_cell_idx];
                    
                        for (uint p2_idx = neighbor_head;            // the neighbors we're calculating
                            p2_idx < neighbor_head + neighbor_count; // forces with
                            p2_idx++)
                        {
                            if (p2_idx == p1_idx) continue; // don't calculate forces with yourself
                            if (p2_idx == next_idx) continue; // don't calculate if your bonded
                            if (p2_idx == prev_idx) continue;

                            // Same filament forces
                            if (p1->filament_id == particles[p2_idx].filament_id)
                            {
                                // normal interaction if two particles in the
                                // same filament are further apart than a min
                                if (abs(int(particles[p2_idx].local_id - p1_local_id)) >= opts.min_local_sep_for_forces)
                                    f += (*particle_particle)(r1, particles[p2_idx].r);
                            }
                            // Different filament forces
                            else
                            {
                                f += (*particle_particle)(r1, particles[p2_idx].r);

                                // TODO: extensile forces
                            }
                        }
                    }

                    p1->f += f;  // apply forces to particle
                }

                // delete force functions
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