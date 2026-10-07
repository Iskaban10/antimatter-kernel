#ifndef ANALYSIS_ENERGY_ANALYSIS_H_
#define ANALYSIS_ENERGY_ANALYSIS_H_

#include "core/real_t.h"

namespace bdm {

class ResourceManager;
class PenningTrapField;

/// Breakdown of total simulation energy, used to validate the Boris
/// integrator (see ComputeTotalEnergy for the physics reasoning).
struct EnergyBreakdown {
  real_t kinetic = 0.0;
  real_t trap_potential = 0.0;
  real_t coulomb_potential = 0.0;
  /// Count of ChargedAntiparticle agents included in this breakdown --
  /// i.e. live, non-annihilated particles at the moment this was
  /// computed. Counted in the same pass as the energy sums (see
  /// ComputeTotalEnergy), not a separate scan.
  int num_active_particles = 0;
  real_t Total() const { return kinetic + trap_potential + coulomb_potential; }
};

/// Computes the total energy of every ChargedAntiparticle currently in
/// `rm`: kinetic (0.5*m*v^2, summed), trap electric potential energy
/// (q*phi(position), summed), and pairwise Coulomb potential energy
/// (k_e*q_i*q_j/r_ij, summed over unique pairs within `coulomb_cutoff_radius`).
///
/// -----------------------------------------------------------------------
/// WHY THIS IS A MEANINGFUL VALIDATION CHECK:
/// -----------------------------------------------------------------------
/// The magnetic field does no work on a charged particle (the magnetic
/// force is always perpendicular to velocity), so it contributes nothing
/// to this energy balance either way -- it's absent from this calculation
/// by construction, not by omission. That means, for a correctly
/// implemented and adequately time-resolved integrator, Total() should
/// stay approximately constant over the course of a simulation. Small
/// oscillation/noise is expected (the Boris pusher is not exactly
/// energy-conserving for the non-magnetic forces, which are applied as a
/// symmetric but still finite-dt kick); a systematic upward or downward
/// DRIFT over many steps is the signal that something is wrong -- e.g.
/// dt too large relative to the cyclotron period (see the dt discussion
/// in params/physics_params.h), or a bug in the Coulomb/Lorentz coupling.
///
/// IMPORTANT: coulomb_cutoff_radius here should be passed the SAME value
/// InteractionParams::cutoff_radius uses for the actual force calculation
/// (interaction/coulomb_interaction.cc). If you pass a different cutoff
/// (e.g. an "infinite" one to get the "true" potential energy), you will
/// see what looks like energy drift whenever a pair of particles crosses
/// the force cutoff's radius -- that would be an artifact of comparing
/// the force model against a different potential model, not a real
/// integration bug. Keeping the two cutoffs consistent is what makes this
/// a clean test of the INTEGRATOR specifically, isolated from the
/// separate, already-documented cutoff-radius approximation itself (see
/// interaction/coulomb_interaction.h design note 1).
///
/// COST: this is an O(N^2) pairwise loop over all currently-alive agents,
/// independent of InteractionParams::cutoff_radius (unlike
/// CoulombInteraction's own neighbor-search-based force calculation).
/// This is intentional -- it's a diagnostic meant to be called
/// periodically (e.g. every N steps), not every step, and a plain pair
/// loop avoids depending on Environment::ForEachNeighbor's exact
/// semantics for a correctness-critical validation tool. At the particle
/// counts this project is currently smoke-testing at (O(10)), the O(N^2)
/// cost is negligible; this will need revisiting if used at much larger N.
EnergyBreakdown ComputeTotalEnergy(ResourceManager* rm,
                                   const PenningTrapField& field,
                                   real_t coulomb_cutoff_radius);

}  // namespace bdm

#endif  // ANALYSIS_ENERGY_ANALYSIS_H_