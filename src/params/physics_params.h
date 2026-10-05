#ifndef PARAMS_PHYSICS_PARAMS_H_
#define PARAMS_PHYSICS_PARAMS_H_

#include <string>

#include "agents/particle_species.h"
#include "core/param/param_group.h"
#include "core/real_t.h"

namespace bdm 
{
  struct PhysicsParams : public ParamGroup 
  {
    BDM_PARAM_GROUP_HEADER(PhysicsParams, 1);

    PhysicsParams(const PhysicsParams&) = default;
    PhysicsParams(PhysicsParams&&) = default;
    PhysicsParams& operator=(const PhysicsParams&) = default;
    PhysicsParams& operator=(PhysicsParams&&) = default;

    // RNG seed for the initial-cloud placement.
    int seed = 42;

    // --- Trap geometry / fields ---------------------------------------------
    // Axial magnetic field strength, Tesla.
    real_t b0 = 1.0;
    // Trap voltage, Volts. See fields/penning_trap.h for sign convention.
    real_t trap_voltage = 10.0;
    // Endcap half-distance, meters.
    real_t z0 = 0.01;
    // Ring electrode radius, meters.
    real_t rho0 = 0.01;

    // --- Time integration ----------------------------------------------------
    // Timestep for the Lorentz push, seconds. Must be small relative to the
    // cyclotron period (2*pi*m / (|q| B0)) -- see forces/lorentz.h.
    real_t dt = 1e-12;
    // Number of simulation steps to run.
    int num_steps = 100000;

    // --- Initial particle cloud ------------------------------------------------
    int num_particles = 10;
    ParticleSpecies species = ParticleSpecies::kPositron;
    real_t initial_cloud_radius = 1e-4;

    // Reads params.json (or `filename`); any key not present keeps the
    // current (default) value. Mirrors CARTopiaX's SimParam::LoadParams.
    void LoadParams(const std::string& filename = "params.json");

    // Prints the effective values -- call after LoadParams to confirm what
    // a given run actually used.
    void PrintParams() const;

    // Throws std::invalid_argument on an obviously invalid configuration.
    void Validate() const;
  };

}

#endif  // PARAMS_PHYSICS_PARAMS_H_