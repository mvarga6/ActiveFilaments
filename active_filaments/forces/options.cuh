#ifndef __AF_FORCES_OTPTIONS_CUH__
#define __AF_FORCES_OTPTIONS_CUH__

#include "../types.h"

namespace af 
{
    enum class FilamentBondType
    {
        None = 0,
        Hookean = 1,
        Fene = 2   
    };

    enum class PairWiseForceType
    {
        None = 0,
        LennardJones = 1,
        WeeksChandlerAnderson = 2
    };

    enum class BondBendingType
    {
        None = 0,
        Cosine = 1,
        Polynomial = 2
    };

    struct ForceKernelOptions
    {
        FilamentBondType backbone_bonds = FilamentBondType::Hookean;
        float backbone_energy_scale = 1.0f;
        float backbone_length_scale = 1.0f;

        PairWiseForceType filament_interaction = PairWiseForceType::LennardJones;
        float interaction_energy_scale = 1.0f;
        float interaction_length_scale = 1.0f;

        BondBendingType filament_bending = BondBendingType::None;
        float bending_energy_scale = 10.0f;
        float bending_length_scale = 1.0f;

        uint min_local_sep_for_forces = 3;
    };
}

#endif