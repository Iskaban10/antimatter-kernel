#ifndef INTERACTION_ANNIHILATION_H_
#define INTERACTION_ANNIHILATION_H_

#include "agents/particle_species.h"
#include "core/behavior/behavior.h"

namespace bdm {

bool IsAntiparticlePair(ParticleSpecies a, ParticleSpecies b);

  class AnnihilationInteraction : public Behavior 
  {
    BDM_BEHAVIOR_HEADER(AnnihilationInteraction, Behavior, 1);

  public:
    AnnihilationInteraction() = default;
    ~AnnihilationInteraction() override = default;

    void Run(Agent* agent) override;
  };

}

#endif