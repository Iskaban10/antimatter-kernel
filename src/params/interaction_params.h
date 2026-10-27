#ifndef PARAMS_INTERACTION_PARAMS_H_
#define PARAMS_INTERACTION_PARAMS_H_

#include <string>

#include "core/param/param_group.h"
#include "core/real_t.h"

namespace bdm 
{
  enum class ComputationMethod 
  {
    kStandard,  // Exact inside the cutoff, zero beyond it.
    kBarnesHut, // Barnes-Hut octree over all particles, does not involve any cut-off radius
    kFastMultipole, // FastMultipole method, to be implemented
  };

  struct EnergyComputationSettings
  {
    ComputationMethod method = ComputationMethod::kStandard;
    real_t cutoff_radius = 0.0;   // for kStandard
    real_t barnes_hut_theta = 0.5;  // for kBarnesHut
    real_t fast_multipole_theta = 0.5;    // for kFastMultipole
    real_t fast_multipole_order = 4;     // for kFastMultipole
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

    std::string coulomb_method = "standard";  // coulomb force computation method
    std::string energy_method = "standard"; // total energy computation method

    real_t cutoff_radius = 0.0; // cut-off radius for the standard calculation method

    real_t barnes_hut_theta = 0.5;  // theta for Barnes-Hut algorithm, lies between 0 and 1. Smaller the value more accurate the analysis

    real_t fast_multipole_order = 4;  // order for fmm, should ideally lie between 3 and 6, higher value leads to more accuracy but more computation time
    real_t fast_multipole_theta = 0.5;  // theta for fast-multipole algorithm, lies between 0 and 1. Smaller the value more accurate the analysis

    ComputationMethod GetCoulombMethod() const;
    ComputationMethod GetEnergyMethod() const;

    void LoadParams(const std::string& filename = "params.json");
    void PrintParams() const;
    void Validate() const;
  };
}
#endif