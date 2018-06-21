#ifndef __AF_NEIGHBOR_FINDER_CUH__
#define __AF_NEIGHBOR_FINDER_CUH__

#include <cuda.h>
#include <thrust/for_each.h>
#include <thrust/transform.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/binary_search.h>
#include <thrust/sequence.h>
#include <thrust/functional.h>
#include <thrust/sort.h>
#include <iostream>

#include "../particle.cuh"
#include "../utilities/vector_type_helpers.cuh"
#include "../utilities/index_map.cuh"
#include "cells.cuh"

namespace af
{
    // Using cell_head_idx and cell_count
    // we can iterate throught all the particles
    // in a particular cell.
    //
    // TODO: Add a cell-by-cell operation to calculate
    // pair-wise forces with all nearby particles for 
    // all particles in that cell
    //  OR 
    // do the same thing particle-by-particle
    thrust::device_vector<uint> particle_cell_ids;
    thrust::device_vector<uint> cell_head_idx;
    thrust::device_vector<uint> cell_tail_idx;
    thrust::device_vector<uint> cell_count;
    bool verbose_neighbor_finding = false;

    class NeighborFinder
    {
        //const float3 size;
        //const uint3 dim;
        Cells cells;
        uint ncells;

    public:

        __host__
        NeighborFinder(Cells cells)
        //: size(cell_size), dim(grid_dim) 
            : cells(cells)
        {
            //ncells = grid_dim.x * grid_dim.y * grid_dim.z;
            this->ncells = cells.count();

            if (cell_head_idx.size() != ncells) 
                cell_head_idx.resize(ncells);
            
            if (cell_tail_idx.size() != ncells) 
                cell_tail_idx.resize(ncells);

            if (cell_count.size() != ncells)
                cell_count.resize(ncells);
        }

        Cells get_cells()
        {
            return cells;
        }

        __host__ __device__
        uint operator()(Particle& p)
        {
            //int3 i = floor_float3(p.r / size);
            //uint X = dim.x;
            //uint XY = X * dim.y;
            //uint cell_id = i.x + i.y*X + i.z*XY;
            //p.cell_id = cell_id; // assign the result to the particle too
            //return cell_id;
            p.cell_id = cells.get_idx(p.r);
            return p.cell_id;
        }

        __host__ __device__
        bool operator()(const Particle& p1, const Particle& p2)
        {
            return (p1.cell_id < p2.cell_id);
        }

        __host__ __device__
        bool operator()(const uint id1, const uint id2)
        {
            return (id1 < id2);
        }

        __host__
        void update(ParticleDeviceArray& particles)
        {
            // match the number of particles
            if (particles.size() != particle_cell_ids.size())
                particle_cell_ids.resize(particles.size());
            
            // get the cell id for each particle
            thrust::transform(particles.begin(), particles.end(), particle_cell_ids.begin(), *this);

            // sort the particle by their cell ids
            thrust::stable_sort_by_key(particle_cell_ids.begin(), particle_cell_ids.end(), particles.begin(), *this);    

            // calculate idx where each cell starts
            thrust::device_vector<uint> cell_ids(ncells);
            thrust::sequence(cell_ids.begin(), cell_ids.end());
            thrust::lower_bound(
                particle_cell_ids.begin(), particle_cell_ids.end(), 
                cell_ids.begin(), cell_ids.end(),
                cell_head_idx.begin()
            );

            // calculate idx where each cell ends
            thrust::upper_bound(
                particle_cell_ids.begin(), particle_cell_ids.end(), 
                cell_ids.begin(), cell_ids.end(),
                cell_tail_idx.begin()
            );

            // count particles in cells (end - begin)
            thrust::transform(
                cell_tail_idx.begin(), cell_tail_idx.end(),
                cell_head_idx.begin(), cell_count.begin(),
                thrust::minus<uint>());

            // recalculate index map
            //ParticleIdxMap idxMap;
            //idxMap.update(particles);

            verbose_print(particles);
        }

    private:
        __host__
        void verbose_print(ParticleDeviceArray& particles)
        {
            // print neighbor list if verbose mode
            if (verbose_neighbor_finding)
            {
                thrust::host_vector<uint> h_cell_heads(cell_head_idx.size());
                thrust::host_vector<uint> h_cell_count(cell_count.size());
                ParticleHostArray h_particles(particles.begin(), particles.end());
                thrust::copy(cell_head_idx.begin(), cell_head_idx.end(), h_cell_heads.begin());
                thrust::copy(cell_count.begin(), cell_count.end(), h_cell_count.begin());

                std::cout << "[cell] pid_1, pid_2, ..." << std::endl;
                for (int i = 0; i < ncells; i++)
                {
                    uint count = h_cell_count[i];
                    if (count > 0)
                    {
                        uint head = h_cell_heads[i];
                        std::cout << "[" << i << "] " << h_particles[head].id;
                        for (int j = head + 1; 
                            j < head + count; 
                            std::cout << ", " << h_particles[j++].id);
                        std::cout << std::endl;
                    }
                }
            }
        }
    };
}

#endif