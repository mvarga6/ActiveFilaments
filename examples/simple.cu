#include "../active_filaments.h"
#include <stdlib.h>
#include <iostream>

using namespace af;

int main(int argc, char *argv[])
{
    // initialize library components
    initialize_aflib(argc, argv);

    // these parameters should be read from cmdline input
    size_t N = 125;
    size_t n = 5;
    size_t n_steps = 10;
    float3 box = make_float3(30,30,0);
    const float dt = 0.0001f;
    const float kBT = 0.1f;
    const float gamma = 0.1f;
    size_t print_rate = 1;

    // temp for testing
    if (argc > 1) N = atoi(argv[1]);
    if (argc > 2) n = atoi(argv[2]);
    if (argc > 3) box.x = atof(argv[3]);
    if (argc > 4) box.y = atof(argv[4]);
    if (argc > 5) box.z = atof(argv[5]);
    if (argc > 6) n_steps = atoi(argv[6]);
    if (argc > 7) print_rate = atoi(argv[7]);
    //if (argc > 7) verbose_neighbor_finding = true;

    // create particle generator from AF library
    //ParticleGeneratorOptions particle_opts;
    //particle_opts.type = InitType::Random;
    //ParticleGenerator particle_generator(particle_opts);

    FilamentGeneratorOptions gen_opts(n);
    gen_opts.bounds = box;
    GridFilamentGenerator filament_generator(gen_opts);

    // generate N particles in a thrust vector
    ParticleHostArray particles = filament_generator.generate(N);

    // build filaments from bonds
    //BondListBuilder bond_builder;
    //bond_builder.add_bond(0,1);
    //bond_builder.add_bond(0,2);
    //bond_builder.add_bond(0,3);
    //bond_builder.add_bond(0,4);

    // setup dynamics objects
    VelocityVerlet velocity_verlet(dt);
    LangevinThermostat langevin_thermostat(gamma, kBT, dt);

    // create things for neighbor finding (will be wrapped by nicer api)
    auto cell_size = make_float3(2,2,2);
    auto grid_dim = make_uint3(box.x/cell_size.x, box.y/cell_size.y, 1);
    Cells cells(grid_dim, cell_size);
    NeighborFinder neighbor_finding(cells);

    // copy them to a gpu array
    ParticleDeviceArray particles_gpu(particles.begin(), particles.end());

    // Force calculations
    ForceKernelOptions force_opts;
    force_opts.filament_interaction = PairWiseForceType::LennardJones;
    force_opts.interaction_energy_scale = 0.25;
    force_opts.interaction_length_scale = 1.0;
    force_opts.backbone_energy_scale = 100.f;
    force_opts.backbone_length_scale = 0.65;
    Forces forces(force_opts, &neighbor_finding);

    // Output Writer
    ParticleXYZWriter output("test", print_rate);

    // Time loop
    std::cout << std::endl << "Starting Simulation Loop..." << std::endl;
    for (int t = 0; t < n_steps; t++)
    {
        //std::cout << t << " ";

        // forces from pair-wise neighbors
        forces.update(particles_gpu);

        // apply langevin noise
        //langevin_thermostat.update(particles_gpu);

        // integrate time
        velocity_verlet.update(particles_gpu);

        // periodically copy back to cpu to print
        //pull_and_print(particles_gpu, particles);

        // write the output files
        ParticleHostArray sample(N*n);
        thrust::copy_n(particles_gpu.begin(), N*n, sample.begin());
        output.update(sample);

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
}