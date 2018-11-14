#ifndef __AF_BODIES_PARTICLE_H__
#define __AF_BODIES_PARTICLE_H__

#include <iostream>
#include <cuda.h>
#include <vector_types.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "utilities.h"

namespace af
{
    class Particle
    {
    public:
        float3 r;     // Position
        float3 v;     // Velocity
        float3 f;     // Current Force Sum
        float3 f_old; // Previous Force Sum

        uint id;          // Unique ID of particle
        uint cell_id;     // ID of the cell particle is in
        uint filament_id; // ID of the filament particle is in
        uint local_id;    // ID of particle within filament
        uint group_id;    // ID of the group of particles

        int next_idx; // idx of next in filament
        int prev_idx; // idx of previous in filament
        float3 t;     // Tangent vector at particle position

        //
        // Class Methods
        //

        __host__
        static Particle create_random()
        {
            Particle p;
            p.r = Random::uniform_float3();
            p.v = zero_float3();
            p.f = zero_float3();
            p.f_old = zero_float3();
            p.id = 0;
            p.cell_id = 0;
            p.filament_id = 0;
            p.group_id = 0;
            p.local_id = 0;
            p.t = zero_float3();
            return p;
        }
    };

    typedef thrust::device_vector<Particle> ParticleDeviceArray;
    typedef thrust::host_vector<Particle> ParticleHostArray;
}

#endif