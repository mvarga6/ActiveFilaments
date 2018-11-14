#ifndef __AF_CELL_SORT_CUH__
#define __AF_CELL_SORT_CUH__

#include <thrust/sort.h>
#include "../particle.cuh"

namespace af
{
    class CellSort
    {
    public:

        __host__ __device__
        CellSort(){}

        __host__ __device__
        bool operator()(const Particle& p1, const Particle& p2)
        {
            return (p1.cell_id < p2.cell_id);
        }

        __host__ __device__
        void update(ParticleDeviceArray& particles)
        {
            thrust::stable_sort(particles.begin(), particles.end(), *this);
        }
    };
}

#endif