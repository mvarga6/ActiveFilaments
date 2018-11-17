#ifndef __AF_NEIGHBOR_FINDER_CUH__
#define __AF_NEIGHBOR_FINDER_CUH__

#include <cuda.h>
#include <thrust/for_each.h>
#include <thrust/transform.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/binary_search.h>
#include <thrust/sequence.h>
#include <thrust/functional.h>
#include <thrust/sort.h>
#include <iostream>

#include "../particle.cuh"
#include "../utilities/vector_type_helpers.cuh"
#include "../utilities/index_map.cuh"
#include "cells.cuh"

#define NOTHEADIDX -1010

namespace af
{
    struct IdxMap
    {
        int from_idx;
        int to_idx;
    };

    struct CreateIdxToHeadMap
    {
        __host__ __device__
        IdxMap operator()(const uint idx, const Particle& p)
        {
            IdxMap _pair;
            _pair.from_idx = idx;
            _pair.to_idx = NOTHEADIDX;
            if (p.local_id == 0)
                _pair.to_idx = p.filament_id;
            return _pair;
        }
    };

    struct KeepHeads
    {
        __host__ __device__
        bool operator()(const IdxMap& map_item)
        {
            return !(map_item.to_idx == NOTHEADIDX);
        }
    };

    struct SelectFromIdx
    {
        __host__ __device__
        uint operator()(const IdxMap& map_item)
        {
            return map_item.from_idx;
        }
    };

    struct SelectToIdx
    {
        __host__ __device__
        uint operator()(const IdxMap& map_item)
        {
            return map_item.to_idx;
        }
    };

    // thrust::device_vector<IdxMap> filament_head_map;

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
    thrust::device_vector<uint> cell_ids;

    thrust::device_vector<IdxMap> filament_idxmap; // destination for particle idx to filament id
    //  |
    //  V
    thrust::device_vector<IdxMap> filament_headonly_idxmap; // destination for only the idx to head filament id
    //  |
    //  V
    thrust::device_vector<uint> filament_head_idx;  // destination for selecting filament head idx
    //  +
    thrust::device_vector<uint> filament_id;        // destination for filament id at that index
    

    bool verbose_neighbor_finding = false;


    __global__ void update_filament_headlist(
        Particle* particles,
        size_t n_particles,
        uint* filament_headlist)
    {
        int pidx = blockIdx.x*blockDim.x + threadIdx.x;
        if (pidx >= n_particles) return // stop if outside of particles
        
        p* = &particles[pidx]
        if (p->local_id != 0) return // stop if not a head

        // Assign the index of this head to the list
        filament_headlist[p->filament_id] = pidx
    }

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

            if (cell_ids.size() != ncells)
                cell_ids.resize(ncells);
        }

        Cells get_cells()
        {
            return cells;
        }

        __host__ __device__
        uint operator()(Particle& p)
        {
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
        void update(ParticleDeviceArray& particles, uint num_filaments)
        {
            this->resize_containers(particles.size(), num_filaments);
            this->sort_particles_by_cells(particles);
            this->count_particles_in_cells();
            this->build_filament_heads_map(particles);

            // recalculate index map
            //ParticleIdxMap idxMap;
            //idxMap.update(particles);

            verbose_print(particles);
        }

    private:

        __host__
        void resize_containers(const uint num_particles, const uint num_filaments)
        {
            // Allocate memory for per particle containers
            if (particle_cell_ids.size() != num_particles)
                particle_cell_ids.resize(num_particles);

            //if (filament_idxmap.size() != num_particles)
            //    filament_idxmap.resize(num_particles);

            // match the number of cells
            // if (ncells != cell_ids.size())
            //     cell_ids.resize(ncells)

            // Allocate memory for per filament containers
            if (filament_head_idx.size() != num_filaments)
                filament_head_idx.resize(num_filaments);

            //if (filament_id.size() != num_filaments)
            //    filament_id.resize(num_filaments);

            //if (filament_headonly_idxmap.size() != num_filaments)
            //    filament_headonly_idxmap.resize(num_filaments);
        }


        __host__
        void sort_particles_by_cells(ParticleDeviceArray& particles)
        {
            // get the cell id for each particle
            thrust::transform(particles.begin(), particles.end(), particle_cell_ids.begin(), *this);

            // sort the particle by their cell ids
            thrust::stable_sort_by_key(particle_cell_ids.begin(), particle_cell_ids.end(), particles.begin(), *this);  
        
            // ParticleHostArray h_particles(particles.begin(), particles.end());
            // for (int i = 0; i < h_particles.size(); i++)
            //     if (h_particles[i].local_id == 0)
            //         std::cout << h_particles[i].filament_id << " @ " << i << std::endl;
        }


        __host__
        void count_particles_in_cells()
        {
            // calculate idx where each cell starts
            // thrust::device_vector<uint> cell_ids(ncells);
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
        }


        __host__
        void build_filament_heads_map(ParticleDeviceArray& particles)
        {
            //thrust::counting_iterator<uint> count_begin(0);
            //thrust::counting_iterator<uint> count_end = count_begin + particles.size();

            // Give the tranformation an enumerated idx so we can create
            // a map from idx to filament id of the only the heads
            //thrust::transform(count_begin, count_end, particles.begin(), filament_idxmap.begin(), CreateIdxToHeadMap());

            // copy the idxmap elements of only the heads
            //auto headonly_end = thrust::remove_copy_if(filament_idxmap.begin(), filament_idxmap.end(), filament_headonly_idxmap.begin(), KeepHeads());

            // Select into two seperate lists so one can sort the other
            //thrust::transform(filament_headonly_idxmap.begin(), headonly_end, filament_head_idx.begin(), SelectFromIdx());
            //thrust::transform(filament_headonly_idxmap.begin(), headonly_end, filament_id.begin(), SelectToIdx());

            // Sort the indices to filament heads BY the filament ids.
            // This will also use to get the idx of a filament head by
            // access that filaments element in the filament_head_idx array.
            //thrust::stable_sort_by_key(filament_id.begin(), filament_id.end(), filament_head_idx.begin());
            const int TPB = 256;
            const int N = particles.size();

            update_filament_headlist<<<TPB,N/TPB + 1>>>
            (
                thrust::raw_pointer_cast(&particles[0]),
                N, thrust::raw_pointer_cast(&filament_head_idx[0])
            );


            // PRINT THE HEADS LIST TO SEE IF IT WORKED!!
            thrust::host_vector<uint> h_fila_heads(filament_head_idx.size());
            thrust::copy(filament_head_idx.begin(), filament_head_idx.end(), h_fila_heads.begin());
            for (int i = 0; i < h_fila_heads.size(); i++)
               std::cout << i << " @ " << h_fila_heads[i] << std::endl;
        }


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