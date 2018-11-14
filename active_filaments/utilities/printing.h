#ifndef __AF_PRINTING_H__
#define __AF_PRINTING_H__

#include <thrust/copy.h>
#include "../particle.cuh"

namespace af
{
    // void pull_and_print(ParticleDeviceArray& from, ParticleHostArray& to)
    // {
    //     std::cout << "Copying & printing data... ";
    //     thrust::copy(from.begin(), from.end(), to.begin());
    //     std::cout << "complete" << std::endl;

    //     for (auto p : to)
    //         std::cout << p.id << " [" << p.cell_id << "] " 
    //                   << p.r.x << " " << p.r.y << " " << p.r.z << " "
    //                   << p.v.x << " " << p.v.y << " " << p.v.z << std::endl;
    // }

    // void pull_and_print(ParticleDeviceArray& from)
    // {
    //     ParticleHostArray to(from.size());
    //     pull_and_print(from, to);
    // }
}

#endif