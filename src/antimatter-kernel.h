// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef ANTIMATTER_KERNEL_H_
#define ANTIMATTER_KERNEL_H_

#include <cmath>
#include <iostream>
#include <random>
#include <chrono>

#include "biodynamo.h"

#include "agents/charged_antiparticle.h"
#include "analysis/energy_analysis.h"
#include "core/param/param.h"
#include "fields/penning_trap.h"
#include "forces/lorentz.h"
#include "interaction/annihilation.h"
#include "interaction/coulomb_interaction.h"
#include "params/interaction_params.h"
#include "params/physics_params.h"
#include "params/stochastic_params.h"

namespace bdm 
{
inline int Simulate(int argc, const char** argv) 
{
  Param::RegisterParamGroup(new PhysicsParams());
  Param::RegisterParamGroup(new InteractionParams());
  Param::RegisterParamGroup(new StochasticParams());

  Simulation simulation(
      argc, argv, [](Param* param) 
      {
        param->Get<PhysicsParams>()->LoadParams("params.json");
        param->Get<InteractionParams>()->LoadParams("params.json");
        param->Get<StochasticParams>()->LoadParams("params.json");
      });
  auto* rm = simulation.GetResourceManager();

  const auto start_time = std::chrono::high_resolution_clock::now();

  const auto* physics = simulation.GetParam()->Get<PhysicsParams>();
  const auto* interaction = simulation.GetParam()->Get<InteractionParams>();
  const auto* stochastic = simulation.GetParam()->Get<StochasticParams>();

  physics->Validate();
  interaction->Validate();
  stochastic->Validate();

  physics->PrintParams();
  interaction->PrintParams();
  stochastic->PrintParams();

  PenningTrapField field = PenningTrapField::MakeFromGeometry( physics->b0, physics->trap_voltage, physics->z0, physics->rho0);

  std::default_random_engine rng(physics->seed);
  std::uniform_real_distribution<real_t> offset( -physics->initial_cloud_radius, physics->initial_cloud_radius);

  for (int i = 0; i < physics->num_particles; ++i) 
  {
    Real3 position{offset(rng), offset(rng), offset(rng)};
    auto* particle = new ChargedAntiparticle(position, physics->species);
    
    if (interaction->coulomb_enabled) { particle->AddBehavior(new CoulombInteraction()); }
    particle->AddBehavior(new LorentzForce(&field));

    if (stochastic->annihilation_enabled) { particle->AddBehavior(new AnnihilationInteraction()); }

    rm->AddAgent(particle);
  }

  const int energy_log_interval = physics->num_steps >= 100 ? physics->num_steps / 100 : 1;

  EnergyBreakdown initial_energy = ComputeTotalEnergy(rm, field, interaction->cutoff_radius);
  std::cout << "\nstep 0  KE=" << initial_energy.kinetic
            << " PE_trap=" << initial_energy.trap_potential
            << " PE_coulomb=" << initial_energy.coulomb_potential
            << " TOTAL=" << initial_energy.Total() << " J" << std::endl;

  for (int step = 1; step <= physics->num_steps; ++step) 
  {
    simulation.GetScheduler()->Simulate(1);

    if (step % energy_log_interval == 0 || step == physics->num_steps) 
    {
      EnergyBreakdown e = ComputeTotalEnergy(rm, field, interaction->cutoff_radius);
      const real_t fractional_drift = initial_energy.Total() != 0.0 ? (e.Total() - initial_energy.Total()) / std::abs(initial_energy.Total()) : 0.0;
      std::cout << "step " << step << "  KE=" << e.kinetic
                << " PE_trap=" << e.trap_potential
                << " PE_coulomb=" << e.coulomb_potential
                << " TOTAL=" << e.Total() << " J"
                << "  (fractional drift from step 0: " << fractional_drift
                << ")" << std::endl;
    }
  }

  rm->ForEachAgent([](Agent* agent) 
  {
    auto* particle = dynamic_cast<ChargedAntiparticle*>(agent);
    if (particle == nullptr) {
      return;
    }
    const auto& pos = particle->GetPosition();
    std::cout << "final position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
  });

  const auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end_time - start_time;
  std::cout << "Simulation completed successfully!" << std::endl;
  std::cout << "Total time taken: " << elapsed.count() << "seconds\n";
  return 0;
}

}
#endif