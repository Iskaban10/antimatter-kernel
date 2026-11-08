#include "forces/lorentz.h"

#include "agents/charged_antiparticle.h"
#include "core/param/param.h"
#include "core/simulation.h"
#include "core/util/math.h"
#include "fields/penning_trap.h"
#include "params/physics_params.h"

namespace bdm
{
  namespace
  {
    // Velocity-only Boris update over a time interval `dt`.
    // v  -->  +half kick  -->  rotate about B  -->  +half kick  -->  v_new
    //         (v_minus)        (v_prime, s_vec)      (v_plus)
    
    Real3 BorisVelocity(const Real3& v, const Real3& force, const Real3& b_field,
                        real_t q, real_t m, real_t dt)
    {
      const real_t half_dt = 0.5 * dt;
      const Real3 half_kick = force * (half_dt / m);

      const Real3 v_minus = v + half_kick;  // intermediate velocity after first-kick

      const Real3 t_vec = b_field * (q * half_dt / m);    // rotation vector in timestep 0.5 * dt
      const Real3 v_prime = v_minus + Math::CrossProduct(v_minus, t_vec);

      const real_t t_mag_sq = t_vec * t_vec;  // dot product
      const Real3 s_vec = t_vec * (2.0 / (1.0 + t_mag_sq)); 

      const Real3 v_plus = v_minus + Math::CrossProduct(v_prime, s_vec);
      return v_plus + half_kick;
    }
  }  // namespace

  void LorentzForce::Run(Agent* agent)
  {
    auto* particle = static_cast<ChargedAntiparticle*>(agent);

    if (particle == nullptr || field_ == nullptr)   return;  // not a charged particle
    if (particle->IsAnnihilated())    return;   // annihilated particle

    const real_t q = particle->GetCharge();
    const real_t m = particle->GetMass();
    if (m <= 0.0 || dt_ <= 0.0)     return;  // misconfigured particle or dt

    const Real3 x = particle->GetPosition();
    Real3 v = particle->GetVelocity();

    const Real3 e_field = field_->GetElectricField(x);
    const Real3 b_field = field_->GetMagneticField(x);
    const Real3 total_force = e_field * q + particle->GetCoulombForce();

    // In the first step only, the stored v is v(0), but the leapfrog expects v(-dt/2).
    if (!velocity_staggered_)
    {
      v = BorisVelocity(v, total_force, b_field, q, m, -0.5 * dt_); // Move it back half a step; the position stays at x(0).
      velocity_staggered_ = true;
    }

    // Normal step: v(n-1/2) -> v(n+1/2), then x(n) -> x(n+1)
    const Real3 v_new = BorisVelocity(v, total_force, b_field, q, m, dt_);
    const Real3 x_new = x + v_new * dt_;

    particle->SetVelocity(v_new);
    particle->SetPosition(x_new);

    // safety fallback if the Coulomb operation skips a step
    particle->SetCoulombForce(Real3{0.0, 0.0, 0.0});
  }
}