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

  // State carried between calls so the potential energy can be time-centered.
  struct EnergyHistory
  {
    bool has_previous = false;
    real_t prev_trap_potential = 0.0;
    real_t prev_coulomb_potential = 0.0;
  };

  // Kinetic energy comes from the stored v(n+1/2) and the potentials from x(n+1).
  // If `history` is given, the potentials are averaged with the previous call's values so they pair with the half-step velocity.
  EnergyBreakdown ComputeTotalEnergy(ResourceManager* rm, const PenningTrapField& field, const EnergyComputationSettings& energy, EnergyHistory* history = nullptr);
}  // namespace bdm

#endif  // ANALYSIS_ENERGY_ANALYSIS_H_