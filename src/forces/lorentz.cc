#include "forces/lorentz.h"

#include "agents/charged_antiparticle.h"
#include "core/param/param.h"
#include "core/simulation.h"
#include "core/util/math.h"
#include "fields/penning_trap.h"
#include "params/physics_params.h"

namespace bdm 
{
  void LorentzForce::Run(Agent* agent) 
  {
    auto* particle = dynamic_cast<ChargedAntiparticle*>(agent);

    if (particle == nullptr || field_ == nullptr) { return; } // no force on if field absent
    if (particle->IsAnnihilated()) { return; }  // annihilated particles are only kept for visualization

    const real_t dt_ = Simulation::GetActive()->GetParam()->Get<PhysicsParams>()->dt;
    const real_t q = particle->GetCharge();
    const real_t m = particle->GetMass();

    if (m <= 0.0 || dt_ <= 0.0) { return; }

    const Real3& x = particle->GetPosition();
    const Real3& v = particle->GetVelocity();

    const Real3 e_field = field_->GetElectricField(x);
    const Real3 b_field = field_->GetMagneticField(x);

    const Real3 coulomb_force = particle->GetCoulombForce();  // defaults to zero in case COulomb force hasn't been run

    // --- Boris pusher, with trap E-field and Coulomb force combined ------
    const real_t half_dt = 0.5 * dt_;
    const Real3 total_force = e_field * q + coulomb_force;
    const Real3 half_kick = total_force * (half_dt / m);

    const Real3 v_minus = v + half_kick;

    const Real3 t_vec = b_field * (q * half_dt / m);
    const Real3 v_prime = v_minus + Math::CrossProduct(v_minus, t_vec);

    const real_t t_mag_sq = t_vec * t_vec;  // Real3 * Real3 is dot product
    const Real3 s_vec = t_vec * (2.0 / (1.0 + t_mag_sq));

    const Real3 v_plus = v_minus + Math::CrossProduct(v_prime, s_vec);
    const Real3 v_new = v_plus + half_kick;

    const Real3 x_new = x + v_new * dt_;

    particle->SetVelocity(v_new);
    particle->SetPosition(x_new);
    particle->SetCoulombForce(Real3{0.0, 0.0, 0.0}); // resets stale force as safety default
  }
}