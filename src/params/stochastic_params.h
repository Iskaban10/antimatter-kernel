#ifndef PARAMS_STOCHASTIC_PARAMS_H_
#define PARAMS_STOCHASTIC_PARAMS_H_

#include <string>

#include "core/param/param_group.h"
#include "core/real_t.h"

namespace bdm 
{
  struct StochasticParams : public ParamGroup 
  {
    BDM_PARAM_GROUP_HEADER(StochasticParams, 1);

    StochasticParams(const StochasticParams&) = default;
    StochasticParams(StochasticParams&&) = default;
    StochasticParams& operator=(const StochasticParams&) = default;
    StochasticParams& operator=(StochasticParams&&) = default;

    bool annihilation_enabled = false;

    real_t contact_radius = 1e-9;
    real_t annihilation_probability = 0.1;

    void LoadParams(const std::string& filename = "params.json");
    void PrintParams() const;
    void Validate() const;
  };
}

#endif