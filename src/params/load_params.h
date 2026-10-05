#ifndef PARAMS_LOAD_PARAMS_H_
#define PARAMS_LOAD_PARAMS_H_

#include <string>

#include "params/interaction_params.h"
#include "params/physics_params.h"
#include "params/stochastic_params.h"

namespace bdm 
{
    // loads physics_params, interaction_params, and stochastic_params from a single JSON config file, 
    void LoadParams(const std::string& path, PhysicsParams* physics, InteractionParams* interaction, StochasticParams* stochastic);
    // print params
    void PrintParams(const PhysicsParams& physics, const InteractionParams& interaction, const StochasticParams& stochastic);
}

#endif