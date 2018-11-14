#ifndef __AF_PAIRWISE_FUNCTOR_FACTORY_CUH__
#define __AF_PAIRWISE_FUNCTOR_FACTORY_CUH__

#include <cuda.h>

#include "pairwise_base.cuh"
#include "../options.cuh"

namespace  af
{
    class PairWiseFactory
    {
    public:
        __host__ __device__
        static PairWiseBase* create(PairWiseForceType type, 
            const float energy_scale, 
            const float length_scale)
        {
            switch(type)
            {
            case PairWiseForceType::LennardJones:
                return new LennardJones(energy_scale, length_scale);
            case PairWiseForceType::WeeksChandlerAnderson:
                return new WeeksChandlerAnderson(energy_scale, length_scale);
            case PairWiseForceType::None:
                return new NoInteraction();
            default:
                return NULL;
            }
        }
    };
}

#endif