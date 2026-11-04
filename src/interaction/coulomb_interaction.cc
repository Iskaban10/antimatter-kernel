#include "interaction/coulomb_interaction.h"

#include <cmath>
#include <cstdint>
#include <vector>

#include "agents/charged_antiparticle.h"
#include "core/container/math_array.h"
#include "core/environment/environment.h"
#include "core/functor.h"
#include "core/resource_manager.h"
#include "core/simulation.h"
#include "utils/constants.h"

namespace bdm 
{
  namespace 
  {
    constexpr real_t kMinSquaredDistance = 1e-30;
  }

  void ComputeStandardCoulombForces(ResourceManager* rm, real_t cutoff_radius) 
  {
    std::vector<ChargedAntiparticle*> particles;
    rm->ForEachAgent([&](Agent* agent) 
    {
      auto* particle = FastCast(agent);
      if (particle != nullptr && !particle->IsAnnihilated())  particles.push_back(particle);
    });

    if (particles.empty())  return;

    if (cutoff_radius <= 0.0) 
    {
      for (auto* particle : particles)  particle->SetCoulombForce(Real3{0.0, 0.0, 0.0});
      return;
    }

    auto* environment = Simulation::GetActive()->GetEnvironment();
    const real_t squared_cutoff = cutoff_radius * cutoff_radius;

    std::vector<Real3> forces(particles.size(), Real3{0.0, 0.0, 0.0});

    const int64_t n = static_cast<int64_t>(particles.size());
    #pragma omp parallel for schedule(dynamic, 64)
    for (int64_t i = 0; i < n; ++i) 
    {
      Agent* agent = particles[i];
      const real_t q_i = particles[i]->GetCharge();
      if (q_i == 0.0)  continue;
      const Real3 pos_i = particles[i]->GetPosition();

      Real3 total_force{0.0, 0.0, 0.0};

      auto accumulate_force = L2F([&](Agent* neighbor, real_t squared_distance) 
      {
        if (neighbor == agent)  return;
        auto* other = FastCast(neighbor);
        if (other == nullptr || other->IsAnnihilated())  return;
        const real_t q_j = other->GetCharge();
        if (q_j == 0.0)  return;
        if (squared_distance < kMinSquaredDistance)  return;  // pair co-incident hence skipped

        const Real3 r_ij = pos_i - other->GetPosition();
        const real_t distance = std::sqrt(squared_distance);
        const real_t force_magnitude = constants::kCoulombConstant * q_i * q_j / squared_distance;
        total_force = total_force + r_ij * (force_magnitude / distance);
      });

      environment->ForEachNeighbor(accumulate_force, *agent, squared_cutoff);
      forces[i] = total_force;
    }

    for (size_t i = 0; i < particles.size(); ++i)  particles[i]->SetCoulombForce(forces[i]);
  }
}