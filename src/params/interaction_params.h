#ifndef PARAMS_INTERACTION_PARAMS_H_
#define PARAMS_INTERACTION_PARAMS_H_

#include <string>

#include "core/param/param_group.h"
#include "core/real_t.h"

namespace bdm 
{
  enum class CoulombMethod 
  {
    kStandard,  // Exact inside the cutoff, zero beyond it.
    kBarnesHut, // Barnes-Hut octree over all particles, does not involve any cut-off radius
    kFastMultipole, // FastMultipole method, to be implemented
  };

  // Coulomb interaction configuration, registered as a ParamGroup
  struct InteractionParams : public ParamGroup 
  {
    BDM_PARAM_GROUP_HEADER(InteractionParams, 1);

    InteractionParams(const InteractionParams&) = default;
    InteractionParams(InteractionParams&&) = default;
    InteractionParams& operator=(const InteractionParams&) = default;
    InteractionParams& operator=(InteractionParams&&) = default;

    bool coulomb_enabled = true;

    std::string coulomb_method = "standard";

    real_t cutoff_radius = 0.0; // cut-off radius for the standard calculation method

    real_t barnes_hut_theta = 0.5;  // theta for Barnes-Hut algorithm, lies between 0 and 1. Smaller the value more accurate the analysis

    CoulombMethod GetCoulombMethod() const;

    void LoadParams(const std::string& filename = "params.json");
    void PrintParams() const;
    void Validate() const;
  };
}
#endif