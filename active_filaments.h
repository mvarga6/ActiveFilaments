#ifndef __ACTIVE_FILAMENTS_H__
#define __ACTIVE_FILAMENTS_H__

#include "active_filaments/particle.cuh"
#include "active_filaments/generators.h"
#include "active_filaments/utilities.h"
#include "active_filaments/integrators.h"
#include "active_filaments/forces.h"
#include "active_filaments/neighbor_finding.h"
#include "active_filaments/io.h"

#include "active_filaments/io/logging.h"
 
namespace af
{
    void initialize_aflib(int argc, char *argv[])
    {
        START_EASYLOGGINGPP(argc, argv);

        LOG(INFO) << "Active Filaments library initialized";
    }
}

#endif