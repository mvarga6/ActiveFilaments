#ifndef __AF_VELOCITY_VERLET_CUH__
#define __AF_vELOCITY_VERLET_CUH__

#include "cuda.h"
#include <thrust/for_each.h>
#include "../particle.cuh"

namespace af
{
    class VelocityVerlet
    {
        const float dt;
        const float dtdt;
    public:

        __host__ __device__
        VelocityVerlet(const float dt)
            : dt(dt), dtdt(dt*dt){}

        __host__ __device__
        void operator()(Particle& p)
        {
            // update position
            p.r += p.v * dt + 0.5f * p.f_old * dtdt;

            // update vecocity
            p.v += (0.5f * (p.f + p.f_old) * dt);

            // save and zero forces
            p.f_old = p.f;
            p.f = zero_float3();
        }

        __host__ __device__
        void update(ParticleDeviceArray& particles)
        {
            thrust::for_each(particles.begin(), particles.end(), *this);
        }

        __host__
        void update(ParticleHostArray& particles)
        {
            thrust::for_each(particles.begin(), particles.end(), *this);
        }
    };
}

#endif