#ifndef ANALYSIS_ENERGY_ANALYSIS_H_
#define ANALYSIS_ENERGY_ANALYSIS_H_

#include "core/real_t.h"
#include <cstddef>
#include "params/interaction_params.h"
#include "interaction/barnes-hut.h"
#include "interaction/fast-multipole.h"

namespace bdm 
{
  class ResourceManager;
  class PenningTrapField;

  // Breakdown of total simulation energy, used to validate the Boris integrator.
  struct EnergyBreakdown 
  {
    real_t kinetic = 0.0;
    real_t trap_potential = 0.0;
    real_t energy_potential = 0.0;
    size_t num_active_particles = 0; // number of non-annihilated particles
    real_t Total() const { return kinetic + trap_potential + energy_potential; }
  };

  EnergyBreakdown ComputeTotalEnergy(ResourceManager* rm, const PenningTrapField& field, const EnergyComputationSettings& energy);
}
#endif