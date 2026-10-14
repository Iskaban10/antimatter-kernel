#include "interaction/coulomb_interaction.h"

#include <cmath>

#include "agents/charged_antiparticle.h"
#include "core/environment/environment.h"
#include "core/functor.h"
#include "core/param/param.h"
#include "core/simulation.h"
#include "core/util/math.h"
#include "params/interaction_params.h"
#include "utils/constants.h"

namespace bdm 
{
  void CoulombInteraction::Run(Agent* agent) 
  {
    auto* particle = FastCast(agent);
    if (particle == nullptr || particle->IsAnnihilated()) { return; }

    const auto* interaction_params = Simulation::GetActive()->GetParam()->Get<InteractionParams>();
    if (!interaction_params->coulomb_enabled) 
    {
      particle->SetCoulombForce(Real3{0.0, 0.0, 0.0});
      return;
    }

    const real_t cutoff_radius = interaction_params->cutoff_radius;
    if (cutoff_radius <= 0.0) 
    {
      particle->SetCoulombForce(Real3{0.0, 0.0, 0.0});
      return;
    }

    const real_t q_i = particle->GetCharge();
    const Real3& pos_i = particle->GetPosition();
    const real_t squared_cutoff = cutoff_radius * cutoff_radius;

    Real3 total_force{0.0, 0.0, 0.0};

    auto accumulate_force = L2F([&](Agent* neighbor, real_t squared_distance) 
    {
      if (neighbor == agent)  return; 
      auto* other = FastCast(neighbor);
      if (other == nullptr || other->IsAnnihilated())  return; 
      const real_t q_j = other->GetCharge();
      if (q_j == 0.0)  return; 

      constexpr real_t kMinSquaredDistance = 1e-30;
      if (squared_distance < kMinSquaredDistance) 
      {
        auto* random = Simulation::GetActive()->GetRandom();
        auto jitter = random->template UniformArray<3>(-1.0, 1.0);
        total_force = total_force + Real3{jitter[0], jitter[1], jitter[2]};
        return;
      }

      const Real3& pos_j = other->GetPosition();
      const Real3 r_ij = pos_i - pos_j;
      const real_t distance = std::sqrt(squared_distance);
      const real_t force_magnitude = constants::kCoulombConstant * q_i * q_j / squared_distance;
      total_force = total_force + r_ij * (force_magnitude / distance);
    });

    auto* environment = Simulation::GetActive()->GetEnvironment();
    environment->ForEachNeighbor(accumulate_force, *agent, squared_cutoff);

    particle->SetCoulombForce(total_force);
  }
}