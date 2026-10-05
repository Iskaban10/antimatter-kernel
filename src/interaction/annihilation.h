#ifndef INTERACTION_ANNIHILATION_H_
#define INTERACTION_ANNIHILATION_H_

#include "agents/particle_species.h"
#include "core/behavior/behavior.h"

namespace bdm {

bool IsAntiparticlePair(ParticleSpecies a, ParticleSpecies b);

/// Simplified stochastic annihilation model. See the conversation turn
/// that introduced this class for full design notes on: (1) the
/// non-physically-sourced flat-probability model, (2) the two-step
/// mark-then-remove pattern matching ChargedAntiparticle::MarkAnnihilated,
/// and (3) the KNOWN, UNRESOLVED thread-safety gap when marking a
/// neighbor's state from inside this agent's Run() under BioDynaMo's
/// (likely parallel) agent execution -- all unchanged by this refactor.
///
/// CHANGED: contact_radius and annihilation_probability are read each step
/// from the registered StochasticParams group rather than being
/// constructor-injected.
class AnnihilationInteraction : public Behavior {
  BDM_BEHAVIOR_HEADER(AnnihilationInteraction, Behavior, 1);

 public:
  AnnihilationInteraction() = default;
  ~AnnihilationInteraction() override = default;

  void Run(Agent* agent) override;
};

}  // namespace bdm

#endif  // INTERACTION_ANNIHILATION_H_