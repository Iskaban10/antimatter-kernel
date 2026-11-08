#include "forces/cooling.h"

#include <cmath>
#include <limits>

#include "agents/charged_antiparticle.h"
#include "utils/constants.h"
#include "core/param/param.h"
#include "core/simulation.h"
#include "core/util/random.h"
#include "fields/penning_trap.h"
#include "params/physics_params.h"

namespace bdm 
{
  real_t CyclotronCooling::CoolingTime(real_t charge, real_t mass, real_t b_mag) 
  {
    const real_t q2 = charge * charge;
    const real_t denom = q2 * q2 * b_mag * b_mag;
    if (denom <= 0.0)  return std::numeric_limits<real_t>::infinity();
    return (0.75 / constants::kCoulombConstant) * mass * mass * mass * constants::kSpeedOfLight * constants::kSpeedOfLight * constants::kSpeedOfLight / denom;
  }

  void CyclotronCooling::Run(Agent* agent) 
  {
    auto* particle = static_cast<ChargedAntiparticle*>(agent);
    if (particle == nullptr || field_ == nullptr)  return;
    if (particle->IsAnnihilated())  return;

    auto* sim = Simulation::GetActive();
    const real_t dt = sim->GetParam()->Get<PhysicsParams>()->dt;

    const real_t q = particle->GetCharge();
    const real_t m = particle->GetMass();
    if (q == 0.0 || m <= 0.0 || dt <= 0.0 || tau_scale_ <= 0.0)  return;

    const Real3 b_field = field_->GetMagneticField(particle->GetPosition());
    const real_t b_mag = std::sqrt(b_field * b_field);
    if (b_mag <= 0.0)  return;

    const real_t tau_E = tau_scale_ * CoolingTime(std::abs(q), m, b_mag);
    if (!std::isfinite(tau_E))  return;

    const Real3 b_hat = b_field * (1.0 / b_mag);
    const Real3 v = particle->GetVelocity();

    const real_t v_par = v * b_hat;           // velocity projection along to magnetic field
    const Real3 v_par_vec = b_hat * v_par;  // conversion back to vector
    Real3 v_perp = v - v_par_vec;           // velocity projection perpendicular to magnetic field

    // Exact exponential damping, stable even when dt >> tau_E.
    const real_t decay = std::exp(-dt / (2.0 * tau_E));
    v_perp = v_perp * decay;

    // Optional thermal floor (fluctuation-dissipation).
    if (bath_temp_K_ > 0.0) 
    {
      const real_t sigma = std::sqrt(constants::kBoltzmannConstant * bath_temp_K_ / m * (1.0 - decay * decay));
      auto* random = sim->GetRandom();
      Real3 noise{random->Gaus(0.0, sigma), random->Gaus(0.0, sigma), random->Gaus(0.0, sigma)};
      noise = noise - b_hat * (noise * b_hat);  // keep noise in the plane perpendicular to B
      v_perp = v_perp + noise;
    }

    particle->SetVelocity(v_par_vec + v_perp);
  }
}