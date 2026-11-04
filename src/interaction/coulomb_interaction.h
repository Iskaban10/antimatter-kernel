#ifndef INTERACTION_COULOMB_INTERACTION_H_
#define INTERACTION_COULOMB_INTERACTION_H_

#include "core/real_t.h"

namespace bdm 
{
  class ResourceManager;
  // Cutoff-based Coulomb force between charged particles, only pairs closer than `cutoff_radius` interact
  void ComputeStandardCoulombForces(ResourceManager* rm, real_t cutoff_radius);
}

#endif