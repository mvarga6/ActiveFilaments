#ifndef __AF_FORCES_PREFORCEKERNEL_CUH__
#define __AF_FORCES_PREFORCEKERNEL_CUH__

#include <cuda.h>
#include <device_launch_parameters.h>
#include <math_functions.h>

#include "../neighbor_finding/cells.cuh"
#include "../particle.cuh"

namespace af 
{
    // launch a thread per cell
    __global__ void preforce_kernel(
        Particle* particles, 
        size_t n_particles,
        uint n_per_filament,
        uint* cell_heads, 
        uint* cell_counts,
        Cells cells,
        uint dimension)
    {
        int cell_idx = blockIdx.x*blockDim.x + threadIdx.x;
        if (cell_idx < cells.count())
        {
            // find particles for this cell
            uint head = cell_heads[cell_idx];
            uint count = cell_counts[cell_idx];

            // how many neighboring cells are there
            uint n_dirs = powf(3, dimension);

            // local variables so were not accessing
            // global memory for every calculation
            float3 r1, r2, t;
            uint p1_local_id;
            float3 next_r, prev_r;
            int next_idx, prev_idx;
            bool is_tail, is_head;

            for (uint p1_idx = head;   // the particles we're applying
                p1_idx < head + count; // forces too
                p1_idx++)
            {
                // get first particle ref
                Particle* p1 = &particles[p1_idx];
                r1 = p1->r;
                p1_local_id = p1->local_id;
                is_tail = p1_local_id == 0;
                is_head = (p1_local_id == n_per_filament - 1);
                next_idx = prev_idx = -1; // TODO: store this value somewhere for access

                // loop over cells
                int searched = 0;
                for (int dir = 0; dir < n_dirs; dir++)
                {
                    searched++;

                    // find particles for the search cell
                    uint neighbor_cell_idx = cells.neighbor_idx(cell_idx, dir);
                    if (neighbor_cell_idx == NO_CELL_IDX) continue;
                    uint cell_head = cell_heads[neighbor_cell_idx];
                    uint cell_count = cell_counts[neighbor_cell_idx];
                
                    for (uint p2_idx = cell_head;        // the neighbor we're calculating
                            p2_idx < cell_head + cell_count; // forces with
                            p2_idx++)
                    {
                        // don't calculate forces with yourself
                        if (p1_idx == p2_idx) continue;
                        
                        // get the ref to other particle
                        Particle* p2 = &particles[p2_idx];
                        r2 = p2->r;
                        
                        // Same filament
                        if (p1->filament_id == p2->filament_id)
                        {
                            // distance along backbone
                            int local_dist = p2->local_id - p1_local_id;

                            if (local_dist == 1)
                            {
                                next_idx = p2_idx;
                                next_r = r2;
                            }
                            else if (local_dist == -1)
                            {
                                prev_idx = p2_idx;
                                prev_r = r2;
                            }
                        }

                        //
                        // Stopping condition check
                        //
                        // only loop neighbors until
                        //  we have the info we need
                        //

                        if ((is_tail && next_idx >= 0) ||
                            (is_head && prev_idx >= 0) ||
                            (next_idx >= 0 && prev_idx >= 0))
                        {
                            // calc tangent
                            if (is_head)
                                t = r1 - prev_r;
                            else if (is_tail)
                                t = next_r - r1;
                            else
                                t = next_r - prev_r;

                            // assign the properties
                            p1->t = norm(t);
                            p1->next_idx = next_idx;
                            p1->prev_idx = prev_idx;

                            // stop the loop on p2
                            p2_idx = cell_head + cell_count;

                            //printf("%d %d %d %d %d %d %d %f\n", dir, neighbor_cell_idx, cell_head, cell_count, searched, next_idx, prev_idx, p1->t);
                        }
                    }
                }
            }
        }
    }
}

#endif