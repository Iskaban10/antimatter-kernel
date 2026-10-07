#ifndef INTERACTION_COULOMB_INTERACTION_H_
#define INTERACTION_COULOMB_INTERACTION_H_

#include "core/behavior/behavior.h"
#include "core/container/math_array.h"

namespace bdm 
{
  class CoulombInteraction : public Behavior 
  {
    BDM_BEHAVIOR_HEADER(CoulombInteraction, Behavior, 1);

  public:
    CoulombInteraction() = default;
    ~CoulombInteraction() override = default;
    void Run(Agent* agent) override;
  };
}

#endif