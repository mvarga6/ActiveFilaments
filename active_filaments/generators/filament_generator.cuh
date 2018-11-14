#ifndef __AF_FILAMENT_GENERATOR_H__
#define __AF_FILAMENT_GENERATOR_H__

#include <cuda.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/generate.h>
#include <vector_functions.h>

#include "../../deps/easylogging++.h"

#include "../particle.cuh"
#include "../defaults.h"

namespace af
{
    struct FilamentGeneratorOptions
    {
        size_t particles_per_filament;
        float bond_length; // initial bond length within filament
        float spacing;  // minimum spacing between filaments
        float3 origin;  // the origin of placement
        float3 bounds;  // relative boundaries to keep particles in 

        FilamentGeneratorOptions(size_t particles_per_filament)
            : particles_per_filament(particles_per_filament),
            bond_length(BOND_LENGTH),
            spacing(SPACING)
            {
                this->origin = make_float3(0,0,0);
                this->bounds = make_float3(1e13,1e13,1e13);
            }
    };



    class FilamentGenerator
    {
    protected:
        uint particle_id; // tracks the particle count
        uint filament_id; // tracks the filament count
        uint local_id; // tracks the id within filament
        float3 R; // current particle position
        FilamentGeneratorOptions opts; // options for generation
        bool full; // space if full and can't place more

    public:
        FilamentGenerator(FilamentGeneratorOptions options)
            : particle_id(0), filament_id(0), local_id(0),
            opts(options), R(options.origin), full(false){}

        virtual Particle operator()() = 0;
        virtual ParticleHostArray generate(size_t n) = 0;
    };



    class GridFilamentGenerator : public FilamentGenerator
    {
    protected:
        uint4 grid;

    public:
        GridFilamentGenerator(FilamentGeneratorOptions options)
            : FilamentGenerator(options){}

        Particle operator()()
        {
            return create_particle();
        }

        ParticleHostArray generate(size_t n_filaments)
        {
            // create host particles with thrust
            size_t np = this->opts.particles_per_filament;
            thrust::host_vector<Particle> parts(n_filaments*np);
            thrust::generate(parts.begin(), parts.end(), *this);
            return parts;
        }

    private:

        Particle create_particle()
        {
            // Build new particle with current state
            Particle p;
            p.r           = this->R;
            p.v           = zero_float3();
            p.f           = zero_float3();
            p.f_old       = zero_float3();
            p.id          = this->particle_id;
            p.filament_id = this->filament_id;
            p.local_id    = this->local_id;

            // TODO: set tangent vector

            this->update();

            return p;
        }

        void increment()
        {
            // increment all values
            this->particle_id = this->particle_id + 1;
            this->filament_id = this->particle_id / opts.particles_per_filament; 
            this->local_id    = this->particle_id % opts.particles_per_filament;
        }

        void update()
        {
            // the values just used for generation
            uint old_filament_id = this->filament_id;
            uint new_particle_id = this->particle_id + 1;
            uint new_filament_id = new_particle_id / opts.particles_per_filament; 
            uint new_local_id    = new_particle_id % opts.particles_per_filament;

            increment();

            // Find placement for next particle
            if (old_filament_id == new_filament_id) 
            {
                // continuing filament along X axis
                R.x += opts.bond_length;
            }
            else // new filament
            {
                // Test if there's room for another filament (in X dir)
                // (NEW HEAD + FILAMENT LENGTH) < XBOUNDS ??
                const float test_tail_rx = (R.x + opts.spacing) + opts.particles_per_filament*opts.bond_length;
                if (test_tail_rx < opts.bounds.x)
                {
                    R.x += opts.spacing; // next is in same row
                    return;
                }

                // Test if there's room for another row (in Y dir)
                const float test_head_ry = (R.y + opts.spacing);
                if (test_head_ry < opts.bounds.y)
                {
                    R.x = opts.origin.x; // back to origin X
                    R.y += opts.spacing; // add a row spacing
                    return;
                }

                const float test_head_rz = R.z + opts.spacing;
                if (test_head_rz < opts.bounds.z)
                {
                    R.x = opts.origin.x; // back to origin X
                    R.y = opts.origin.y; // back to origin Y
                    R.z += opts.spacing; // add depth spacing
                    return;
                }

                LOG(ERROR) << "Simulation space is completely full.";
                LOG(INFO) << "Please retry with less filaments or a large box.";
                LOG(INFO) << "particle: " << new_particle_id;
                LOG(INFO) << "filament: " << new_filament_id;
                LOG(INFO) << "local:    " << new_local_id;
                exit(1);
            }

            if (!(R <= opts.bounds).all())
            {
                LOG(ERROR) << "Placed particle outside of simulation box.";
                LOG(INFO) << "Please retry with less filaments or a large box.";
                LOG(INFO) << "particle: " << new_particle_id;
                LOG(INFO) << "filament: " << new_filament_id;
                LOG(INFO) << "local:    " << new_local_id;
                exit(1);
            }
        }
    };
}

#endif