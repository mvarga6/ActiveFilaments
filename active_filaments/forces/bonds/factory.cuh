#ifndef __AF_BOND_FUNCTOR_FACTORY_CUH__
#define __AF_BOND_FUNCTOR_FACTORY_CUH__

#include <cuda.h>

#include "hookean.cuh"
#include "fene.cuh"
#include "../options.cuh"

namespace af
{
    class BondFactory
    {
    public:
        __host__ __device__
        static BondBase* create(FilamentBondType type, 
            const float energy_scale, 
            const float length_scale)
        {
            switch(type)
            {
            case FilamentBondType::None:
                return new NoBond();
            case FilamentBondType::Hookean:
                return new HookeanBond(energy_scale, length_scale);
            case FilamentBondType::Fene:
                return new FeneBond(energy_scale, length_scale);
            default:
                return NULL;
            }
        }

        __host__ __device__ 
        static BendingBase* create(BondBendingType type,
            const float bending_energy)
        {
            switch(type)
            {
            case BondBendingType::None:
                return new CosineBending(0);
            case BondBendingType::Cosine:
                return new CosineBending(bending_energy);
            default:
                return NULL;
            }
        }
    };
}

#endif