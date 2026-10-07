#ifndef PARAMS_INTERACTION_PARAMS_H_
#define PARAMS_INTERACTION_PARAMS_H_

#include <string>

#include "core/param/param_group.h"
#include "core/real_t.h"

namespace bdm 
{
  // Coulomb interaction configuration, registered as a ParamGroup
  struct InteractionParams : public ParamGroup 
  {
    BDM_PARAM_GROUP_HEADER(InteractionParams, 1);

    InteractionParams(const InteractionParams&) = default;
    InteractionParams(InteractionParams&&) = default;
    InteractionParams& operator=(const InteractionParams&) = default;
    InteractionParams& operator=(InteractionParams&&) = default;

    bool coulomb_enabled = true;
    real_t cutoff_radius = 0.0;

    void LoadParams(const std::string& filename = "params.json");
    void PrintParams() const;
    void Validate() const;
  };
}
#endif