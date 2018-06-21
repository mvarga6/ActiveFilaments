#ifndef __AF_FORCES_NOISE_H__
#define __AF_FORCES_NOISE_H__

#include <cuda.h>
#include <thrust/for_each.h>
#include <vector_functions.h>

#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/normal_distribution.h>

#include "../particle.cuh"
#include "../utilities/random.cuh"

namespace af
{
    class LangevinThermostat
    {
        float noise;
        float g;

        thrust::minstd_rand rng;
        thrust::random::normal_distribution<float> dist;

    public:

        __host__ __device__
        LangevinThermostat(const float gamma, const float kBT, const float dt)
            : g(gamma)
        {
            this->noise = sqrt(2.0f * gamma * kBT / dt);
        }

        __host__ __device__
        void operator()(Particle& p)
        {
            // random term TODO: get RNG to work
            //float Rx = dist(rng);
            //float Ry = dist(rng);
            //float Rz = dist(rng);
            //p.f += (noise * make_float3(Rx, Ry, Rz, 0));

            // drag term
            p.f -= (g * p.v);
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