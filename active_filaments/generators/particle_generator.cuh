#ifndef __AF_PARTICLE_GENERATOR_H__
#define __AF_PARTICLE_GENERATOR_H__

#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/generate.h>
#include "../particle.cuh"

namespace af
{
    enum class InitType
    {
        Random,
        Linear,
        Custom
    };

    struct ParticleGeneratorOptions
    {
        InitType type;
        uint filament_length;
    };

    class ParticleGenerator
    {
        ParticleGeneratorOptions opts;
        uint particle_id; // tracks the particle count
        uint filament_id; // tracks the filament count
        uint local_id; // tracks the id within filament
    public:
        ParticleGenerator(ParticleGeneratorOptions options)
            : opts(options), particle_id(0),
            filament_id(0), local_id(0){}

        Particle operator()()
        {
            Particle p;

            // construction
            switch (opts.type)
            {
            case InitType::Random:
                p = Particle::create_random();
                break;

            case InitType::Linear:
                
                break;

            case InitType::Custom:

                break;
            }

            p.id = particle_id++;
            return p;
        }

        ParticleHostArray generate(size_t n_particles)
        {
            // create host particles with thrust
            thrust::host_vector<Particle> parts(n_particles);
            thrust::generate(parts.begin(), parts.end(), *this);
            return parts;
        }


    };
}

#endif