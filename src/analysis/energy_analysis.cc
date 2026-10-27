#include "analysis/energy_analysis.h"

#include <cmath>
#include <cstdint>
#include <vector>

#include "agents/charged_antiparticle.h"
#include "core/resource_manager.h"
#include "fields/penning_trap.h"
#include "utils/constants.h"

namespace bdm
{
  namespace
  {
    real_t ComputeStandardEnergy(const std::vector<ChargedAntiparticle*>& particles, real_t cutoff)
    {
      const real_t squared_cutoff = cutoff * cutoff;
      const int64_t n = static_cast<int64_t>(particles.size());
      real_t sum = 0.0;

      #pragma omp parallel for reduction(+ : sum) schedule(dynamic, 16)
      for (int64_t i = 0; i < n; ++i)
      {
        const real_t q_i = particles[i]->GetCharge();
        if (q_i == 0.0) continue;  // neutral species contribute no Coulomb PE
        const Real3 x_i = particles[i]->GetPosition();

        for (int64_t j = i + 1; j < n; ++j)
        {
          const real_t q_j = particles[j]->GetCharge();
          if (q_j == 0.0) continue;

          const Real3 r_ij = x_i - particles[j]->GetPosition();
          const real_t d2 = r_ij * r_ij;
          if (d2 > squared_cutoff) continue;  // outside the force model's cutoff
          if (d2 < 1e-30) continue;           // degenerate coincident pair

          sum += constants::kCoulombConstant * q_i * q_j * (1.0 / std::sqrt(d2) - 1.0 / cutoff);
        }
      }
      return sum;
    }

    real_t ComputeBarnesHutEnergy(const std::vector<ChargedAntiparticle*>& particles, real_t theta)
    {
      const size_t n = particles.size();
      std::vector<BarnesHutTree::Source> sources(n);
      for (size_t i = 0; i < n; ++i)
      {
        sources[i].position = particles[i]->GetPosition();
        sources[i].charge = particles[i]->GetCharge();
      }

      BarnesHutTree tree(theta);
      tree.Build(sources);

      real_t sum = 0.0;
      const int64_t count = static_cast<int64_t>(n);

      #pragma omp parallel for reduction(+ : sum) schedule(static)
      for (int64_t i = 0; i < count; ++i)
      {
        const size_t idx = static_cast<size_t>(i);
        if (sources[idx].charge == 0.0) continue;
        sum += sources[idx].charge * tree.ComputePotential(sources[idx].position, idx);
      }
      return 0.5 * constants::kCoulombConstant * sum; // mutiply by half since each pair is taken twice
    }

    real_t ComputeFastMultipoleEnergy(const std::vector<ChargedAntiparticle*>& particles, int order, real_t theta)
    {
      std::vector<FastMultipoleSolver::Source> sources(particles.size());
      for (size_t i = 0; i < particles.size(); ++i)
      {
        sources[i].position = particles[i]->GetPosition();
        sources[i].charge = particles[i]->GetCharge();
      }
      FastMultipoleSolver solver(order, theta);
      return solver.ComputePotentialEnergy(sources);
    }
  }  // namespace

  EnergyBreakdown ComputeTotalEnergy(ResourceManager* rm, const PenningTrapField& field, const EnergyComputationSettings& energy)
  {
    EnergyBreakdown result;

    // Collect non-annihilated particles so the solver can index directly since annihilated particle is not part of system's energy.
    std::vector<ChargedAntiparticle*> particles;
    rm->ForEachAgent([&](Agent* agent)
    {
      auto* particle = static_cast<ChargedAntiparticle*>(agent);
      if (particle != nullptr && !particle->IsAnnihilated()) particles.push_back(particle);
    });

    result.num_active_particles = particles.size();

    for (auto* particle : particles)
    {
      const real_t m = particle->GetMass();
      const real_t q = particle->GetCharge();
      const Real3& v = particle->GetVelocity();
      const Real3& x = particle->GetPosition();

      result.kinetic += 0.5 * m * (v * v);
      result.trap_potential += q * field.GetPotential(x);
    }

    // Coulomb potential energy: solver chosen by the user.
    switch (energy.method)
    {
      case ComputationMethod::kFastMultipole:
        result.energy_potential = ComputeFastMultipoleEnergy(particles, energy.fast_multipole_order, energy.fast_multipole_theta);
        break;
      case ComputationMethod::kBarnesHut:
        result.energy_potential = ComputeBarnesHutEnergy(particles, energy.barnes_hut_theta);
        break;
      case ComputationMethod::kStandard:
      default:
        result.energy_potential = ComputeStandardEnergy(particles, energy.cutoff_radius);
        break;
    }

    return result;
  }
}  // namespace bdm