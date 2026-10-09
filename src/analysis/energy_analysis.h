#ifndef ANALYSIS_ENERGY_ANALYSIS_H_
#define ANALYSIS_ENERGY_ANALYSIS_H_

#include "core/real_t.h"

namespace bdm 
{
  class ResourceManager;
  class PenningTrapField;

  // Breakdown of total simulation energy, used to validate the Boris integrator.
  struct EnergyBreakdown 
  {
    real_t kinetic = 0.0;
    real_t trap_potential = 0.0;
    real_t coulomb_potential = 0.0;
    int num_active_particles = 0; // no. of non-annihilated particles
    real_t Total() const { return kinetic + trap_potential + coulomb_potential; }
  };

  EnergyBreakdown ComputeTotalEnergy(ResourceManager* rm, const PenningTrapField& field, real_t coulomb_cutoff_radius);
}
#endif