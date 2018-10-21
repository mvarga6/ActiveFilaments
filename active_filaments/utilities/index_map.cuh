#ifndef __AF_UTILITIES_INDEX_MAP_CUH__
#define __AF_UTILITIES_INDEX_MAP_CUH__

#include <cuda_runtime.h>
#include "../particle.cuh"
#include "enumorate.cuh"

namespace af
{
    thrust::device_vector<uint> id_to_idx;

    class ParticleIdxMap
    {
    public:

        __host__ __device__
        ParticleIdxMap(){}

        __host__ __device__
        void operator()(const enumeratedParticle& id_particle)
        {
            // the current index of this particle
            uint idx = id_particle.get<0>();

            // get the particle's id
            uint pid = id_particle.get<1>().id;

            // map the particle id to it's index
            //id_to_idx[pid] = idx;
        }

        __host__
        void update(ParticleDeviceArray& particles)
        {
            EnumeratedDeviceParticles parts(particles);
            thrust::for_each(parts.begin(), parts.end(), *this);
        }

        // __host__
        // static void update(ParticleDeviceArray& particles)
        // {
        //     EnumoratedDeviceParticles parts(particles);
        //     thrust::for_each(parts.begin(), parts.end(), ParticleIdxMap());
        // }
    };
}

#endif