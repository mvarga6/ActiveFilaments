#include "../active_filaments.h"
#include <stdlib.h>
#include <iostream>

#include "../deps/easylogging++.h"

using namespace af;

int main(int argc, char *argv[])
{
    // initialize library components
    initialize_aflib(argc, argv);

    // these parameters should be read from cmdline input
    size_t num_filaments = 125;
    size_t particles_per_filament = 5;
    size_t num_steps = 100;
    float3 box = make_float3(30,30,0);
    const float dt = 0.0001f;
    const float kBT = 0.1f;
    const float gamma = 0.1f;
    size_t print_rate = 1;

    // temp for testing
    if (argc > 1) num_filaments = atoi(argv[1]);
    if (argc > 2) particles_per_filament = atoi(argv[2]);
    if (argc > 3) box.x = atof(argv[3]);
    if (argc > 4) box.y = atof(argv[4]);
    if (argc > 5) box.z = atof(argv[5]);
    if (argc > 6) num_steps = atoi(argv[6]);
    if (argc > 7) print_rate = atoi(argv[7]);
    if (argc > 8) verbose_neighbor_finding = bool(atoi(argv[8]));
    LOG(INFO) << "Loaded parameters.";

    // create particle generator from AF library
    //ParticleGeneratorOptions particle_opts;
    //particle_opts.type = InitType::Random;
    //ParticleGenerator particle_generator(particle_opts);


    FilamentGeneratorOptions gen_opts(particles_per_filament);
    gen_opts.bounds = box;
    gen_opts.spacing = 0.8;
    gen_opts.bond_length = 0.9;
    GridFilamentGenerator filament_generator(gen_opts);
    LOG(INFO) << "FilamentGenerator created.";

    // generate N filaments in a thrust vector
    ParticleHostArray particles = filament_generator.generate(num_filaments);
    LOG(INFO) << "ParticleHostArray created.";

    // build filaments from bonds
    //BondListBuilder bond_builder;
    //bond_builder.add_bond(0,1);
    //bond_builder.add_bond(0,2);
    //bond_builder.add_bond(0,3);
    //bond_builder.add_bond(0,4);

    // setup dynamics objects
    VelocityVerlet velocity_verlet(dt);
    LOG(INFO) << "VelocityVerlet integrator created.";

    LangevinThermostat langevin_thermostat(gamma, kBT, dt);
    LOG(INFO) << "LangevinThermostat created.";

    // create things for neighbor finding (will be wrapped by nicer api)
    auto cell_size = make_float3(1.5,1.5,1.5);
    auto grid_dim = make_uint3(box.x/cell_size.x, box.y/cell_size.y, 1);
    Cells cells(grid_dim, cell_size);
    NeighborFinder neighbor_finding(cells);
    LOG(INFO) << "NeighborFinder created.";

    // copy them to a gpu array
    ParticleDeviceArray particles_gpu(particles.begin(), particles.end());
    LOG(INFO) << "ParticleDeviceArray created.";

    // Force calculations
    ForceKernelOptions force_opts;
    force_opts.filament_interaction = PairWiseForceType::LennardJones;
    force_opts.interaction_energy_scale = 0.4;
    force_opts.interaction_length_scale = 1.0;
    force_opts.backbone_energy_scale = 100.f;
    force_opts.backbone_length_scale = 0.65;
    Forces forces(force_opts, &neighbor_finding);
    LOG(INFO) << "ForceKernel created.";

    // Output Writer
    ParticleXYZWriter output("test", print_rate);
    LOG(INFO) << "ParticleWriter created.";

    // Time loop
    std::cout << std::endl << "Starting Simulation Loop..." << std::endl;
    for (int t = 0; t < num_steps; t++)
    {
        try
        {
            //std::cout << t << " ";

            // forces from pair-wise neighbors
            forces.update(particles_gpu, particles_per_filament, num_filaments); //TODO: pass n in better way

            // apply langevin noise
            //langevin_thermostat.update(particles_gpu);

            // integrate time
            velocity_verlet.update(particles_gpu);

            // write the output files
            output.update(particles_gpu);

            // print first 5 particles
            if (verbose_neighbor_finding)
            {
                ParticleHostArray sample(5);
                thrust::copy_n(particles_gpu.begin(), 5, sample.begin());

                for (auto p : sample)
                    std::cout << p.id << " [" << p.cell_id << "] " 
                        << p.r.x << " " << p.r.y << " " << p.r.z << " "
                        << p.v.x << " " << p.v.y << " " << p.v.z << std::endl;
            }
        }
        catch(...)
        {
            LOG(ERROR) << "Fatal error encountered in simulation. :(";
            break;
        }
        
    }
}