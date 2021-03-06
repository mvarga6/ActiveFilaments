#ifndef __AF_CELLS_CUH__
#define __AF_CELLS_CUH__

#include <cuda.h>
#include <cuda_runtime.h>
#include "../utilities/vector_type_helpers.cuh"

#include "../types.h"

#define NO_CELL_IDX 654321

namespace af 
{
    struct Cells
    {
        uint3 grid;
        float3 size;
        bool pbc;
    
        __host__ __device__
        Cells(uint3 cell_grid, float3 cell_size, bool periodic = true)
            : grid(cell_grid), size(cell_size), pbc(periodic){}

        __host__ __device__
        uint count()
        {
            return grid.x * grid.y * grid.z;
        }

        __host__ __device__
        uint neighbor_idx(uint idx, uint dir)
        {
            uint3 ijk = get_ijk(idx);
            int3 nijk = ijk + dijk(dir);

            if (pbc) apply_pbc(nijk);

            if (!cell_exists(nijk))
                return NO_CELL_IDX;

            return get_idx(nijk);
        }

        __host__ __device__
        int3 dijk(uint dir)
        {
            int di = unit_cycle(dir % 3);
            int dj = unit_cycle((dir / 3) % 3);
            int dk = unit_cycle(dir / 9);
            return make_int3(di, dj, dk);
        }

        // maps [0, 1, 2] -> [0, 1, -1]
        __host__ __device__
        int unit_cycle(uint i)
        {
            if (i > 1)
            {
                return -1;
            }
            else return i;
        }

        __host__ __device__
        uint3 get_ijk(uint idx)
        {
            uint X = grid.x;
            uint XY = grid.y * X;
            return make_uint3(
                idx % X,
                idx / X,
                idx / XY
            );
        }

        __host__ __device__
        uint get_idx(uint3 ijk)
        {
            uint X = grid.x;
            uint XY = grid.y * X;
            return ijk.x + ijk.y*X + ijk.z*XY;
        }

        __host__ __device__
        uint get_idx(int3 ijk)
        {
            uint X = grid.x;
            uint XY = grid.y * X;
            return ijk.x + ijk.y*X + ijk.z*XY;
        }

        __host__ __device__
        uint get_idx(float3 pos)
        {
            return get_idx(floor_float3(pos / size));
        }

        __host__ __device__
        void apply_pbc(int3& ijk)
        {
            ijk = (grid + ijk) % grid;
        }

        __host__ __device__
        bool cell_exists(int3& ijk)
        {
            return in_range(ijk, zero_uint3(), grid).all();
        }
    };
}

#endif