#include "analysis/energy_analysis.h"

#include <cmath>
#include <vector>

#include "agents/charged_antiparticle.h"
#include "core/resource_manager.h"
#include "fields/penning_trap.h"
#include "utils/constants.h"

namespace bdm {

EnergyBreakdown ComputeTotalEnergy(ResourceManager* rm,
                                   const PenningTrapField& field,
                                   real_t coulomb_cutoff_radius) {
  EnergyBreakdown result;

  // Collect live, non-annihilated particles first (rather than computing
  // inline inside ForEachAgent) so the O(N^2) Coulomb pair loop below can
  // index directly, and so annihilated agents (still present for one step
  // per ChargedAntiparticle::MarkAnnihilated's documented contract) are
  // cleanly excluded from the energy balance -- an annihilated particle's
  // kinetic/potential energy is not meaningfully part of "the system's
  // energy" for this check.
  std::vector<ChargedAntiparticle*> particles;
  rm->ForEachAgent([&](Agent* agent) {
    auto* particle = dynamic_cast<ChargedAntiparticle*>(agent);
    if (particle != nullptr && !particle->IsAnnihilated()) {
      particles.push_back(particle);
    }
  });

  result.num_active_particles = static_cast<int>(particles.size());

  for (auto* particle : particles) {
    const real_t m = particle->GetMass();
    const real_t q = particle->GetCharge();
    const Real3& v = particle->GetVelocity();
    const Real3& x = particle->GetPosition();

    const real_t speed_squared = v * v;  // Real3 * Real3 is dot product
    result.kinetic += 0.5 * m * speed_squared;
    result.trap_potential += q * field.GetPotential(x);
  }

  const real_t squared_cutoff = coulomb_cutoff_radius * coulomb_cutoff_radius;
  for (size_t i = 0; i < particles.size(); ++i) {
    for (size_t j = i + 1; j < particles.size(); ++j) {
      const real_t q_i = particles[i]->GetCharge();
      const real_t q_j = particles[j]->GetCharge();
      if (q_i == 0.0 || q_j == 0.0) {
        continue;  // neutral species contribute no Coulomb PE
      }

      const Real3 r_ij = particles[i]->GetPosition() - particles[j]->GetPosition();
      const real_t squared_distance = r_ij * r_ij;
      if (squared_distance > squared_cutoff) {
        continue;  // outside the force model's cutoff -- see header note
      }
      if (squared_distance < 1e-30) {
        continue;  // degenerate coincident pair; skip rather than divide by ~0
      }

      const real_t distance = std::sqrt(squared_distance);
      result.coulomb_potential += constants::kCoulombConstant * q_i * q_j / distance;
    }
  }

  return result;
}

}