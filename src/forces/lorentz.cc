#include "forces/lorentz.h"
#include "agents/charged_antiparticle.h"
#include "fields/penning_trap.h"
#include "params/physics_params.h"
#include "core/simulation.h"
#include "core/param/param.h"

namespace bdm {

BDM_BEHAVIOR_HEADER_IMPL(LorentzForce);

void LorentzForce::Run(Agent* agent) 
{
  auto* particle = dynamic_cast<ChargedAntiparticle*>(agent);
  if (!particle) return;

  auto* sim = Simulation::GetActive();
  auto* pparams = sim->GetParam()->Get<PhysicsParams>();

  PenningTrap trap(pparams->trap_v0, pparams->trap_r0, pparams->trap_z0, pparams->b_field_z);

  double dt = sim->GetParam()->simulation_time_step;
  double q = particle->GetCharge();
  double m = particle->GetMass();

  Real3 pos = particle->GetPosition();
  Real3 vel = particle->GetVelocity();

  Real3 E = trap.GetElectricField(pos);
  Real3 B = trap.GetMagneticField(pos);

  // Lorentz Force: F = q * (E + v x B)
  Real3 v_cross_B = {
      vel[1] * B[2] - vel[2] * B[1],
      vel[2] * B[0] - vel[0] * B[2],
      vel[0] * B[1] - vel[1] * B[0]
  };

  Real3 force = {
      q * (E[0] + v_cross_B[0]),
      q * (E[1] + v_cross_B[1]),
      q * (E[2] + v_cross_B[2])
  };

  Real3 accel = {force[0] / m, force[1] / m, force[2] / m}; // accn = force / mass

  // Update velocity and position
  particle->SetVelocity({vel[0] + accel[0] * dt, vel[1] + accel[1] * dt, vel[2] + accel[2] * dt});
  particle->SetPosition({pos[0] + vel[0] * dt, pos[1] + vel[1] * dt, pos[2] + vel[2] * dt});
}

} // namespace bdm