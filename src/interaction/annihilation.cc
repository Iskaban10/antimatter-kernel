#include "interaction/annihilation.h"

#include "agents/charged_antiparticle.h"
#include "core/environment/environment.h"
#include "core/functor.h"
#include "core/param/param.h"
#include "core/simulation.h"
#include "params/stochastic_params.h"

namespace bdm 
{
  bool IsAntiparticlePair(ParticleSpecies a, ParticleSpecies b) 
  {
    auto is_pair = [](ParticleSpecies x, ParticleSpecies y, ParticleSpecies p, ParticleSpecies q) 
    {
      return (x == p && y == q) || (x == q && y == p);
    };
    return is_pair(a, b, ParticleSpecies::kAntiproton, ParticleSpecies::kProton) || is_pair(a, b, ParticleSpecies::kPositron, ParticleSpecies::kElectron);
  }

  void AnnihilationInteraction::Run(Agent* agent) 
  {
    auto* particle = FastCast(agent);
    if (particle == nullptr) { return; }

    const auto* stochastic_params = Simulation::GetActive()->GetParam()->Get<StochasticParams>();

    if (!stochastic_params->annihilation_enabled) { return; } // return if annihilation not enabled
    if (particle->IsAnnihilated()) 
    {
      particle->RemoveFromSimulation();
      return;
    }

    const real_t contact_radius = stochastic_params->contact_radius;
    if (contact_radius <= 0.0) { return; }
    const real_t annihilation_probability = stochastic_params->annihilation_probability;

    const ParticleSpecies species_i = particle->GetSpecies();
    const real_t squared_contact = contact_radius * contact_radius;
    bool paired_this_step = false;

    auto check_partner = L2F([&](Agent* neighbor, real_t squared_distance) 
    {
      if (paired_this_step || neighbor == agent) { return; }
      auto* other = FastCast(neighbor);
      if (other == nullptr || other->IsAnnihilated()) { return; }
      if (!IsAntiparticlePair(species_i, other->GetSpecies())) { return; }

      auto* random = Simulation::GetActive()->GetRandom();
      if (random->Uniform(0.0, 1.0) < annihilation_probability) 
      {
        particle->MarkAnnihilated();
        other->MarkAnnihilated();
        paired_this_step = true;
      }
    });

    auto* environment = Simulation::GetActive()->GetEnvironment();
    environment->ForEachNeighbor(check_partner, *agent, squared_contact);
  }
}