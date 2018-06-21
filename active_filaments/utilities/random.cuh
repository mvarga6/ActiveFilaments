#ifndef __AF_RANDOM_CUH__
#define __AF_RANDOM_CUH__

#include <vector_types.h>
#include <vector_functions.hpp>
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/normal_distribution.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>

namespace af
{
    class Random
    {
    public:

        static bool initialized;

        static void initialize()
        {
            std::cout << "Initializing Random module" << std::endl;
            srand(time(NULL));
            initialized = true;
        }

        static float uniform_float(const float max = 1.0f, const float min = 0.0f)
        {
            if (!Random::initialized) initialize();
            return float(rand() % 10000) / 10000.f;
        }

        static float3 uniform_float3(const float max = 1.0f, const float min = 0.0f)
        {
            return make_float3(
                uniform_float(max, min), 
                uniform_float(max, min),
                uniform_float(max, min)
            );
        }

        __host__ __device__
        static float3 gaussian_float3()
        {
            // This should be stored somewhere else
            thrust::minstd_rand rng;
            thrust::random::normal_distribution<float> dist;
            
            float Rx = dist(rng);
            float Ry = dist(rng);
            float Rz = dist(rng);
            return make_float3(Rx, Ry, Rz);
        }
    };

    bool Random::initialized = false;
}

#endif