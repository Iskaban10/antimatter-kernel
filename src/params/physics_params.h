#ifndef PARAMS_PHYSICS_PARAMS_H_
#define PARAMS_PHYSICS_PARAMS_H_

#include <map>
#include <vector>

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

    // Timestep for the Lorentz push, seconds, should be small relative to the cyclotron period (2*pi*m / (|q| B0))
    real_t dt = 1e-9;
    size_t num_steps = 5000; // Number of simulation steps to run.

    // Cooling Params
    real_t bath_temp = 0.5;
    real_t cooling_scale = 10;

    // --- Initial particle cloud ------------------------------------------------
    int num_particles = 10;
    std::map<ParticleSpecies, size_t> populations = {{ParticleSpecies::kPositron, 10}};
    std::map<ParticleSpecies, real_t> particle_temperatures = {{ParticleSpecies::kElectron, 4},{ParticleSpecies::kPositron, 4},{ParticleSpecies::kAntiproton, 1000}}; // initial temperatures of particles
    Real3 initial_cloud_radii{0.0001, 0.0001, 0.0001};  // radii along each axis for the cloud

    void LoadParams(const std::string& filename = "params.json"); // loads params.json
    void PrintParams() const; // print params that have been loaded

    // Throws std::invalid_argument on an obviously invalid configuration.
    void Validate() const;
  };
}
#endif