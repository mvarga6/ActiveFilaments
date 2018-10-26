#ifndef __AF_CELL_ASSIGN_CUH__
#define __AF_CELL_ASSIGN_CUH__

#include "../particle.cuh"
#include "../utilities/vector_type_helpers.cuh"

namespace af
{
    class CellAssign
    {
        const float3 size;
        const uint3 dim;
    public:

        __host__ __device__
        CellAssign(const float3 cell_size, const uint3 grid_dim)
            : size(cell_size), dim(grid_dim){}

        __host__ __device__
        uint operator()(Particle& p)
        {
            int3 i = floor_float3(p.r / size);
            unsigned int X = dim.x;
            unsigned int XY = X * dim.y;
            return i.x + i.y*X + i.z*XY;
        }
    };
}

#endif