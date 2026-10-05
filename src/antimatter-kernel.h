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

#include <iostream>
#include <random>

#include "biodynamo.h"

#include "agents/charged_antiparticle.h"
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

    Simulation simulation( argc, argv, [](Param* param) 
        {
          param->Get<PhysicsParams>()->LoadParams("params.json");
          param->Get<InteractionParams>()->LoadParams("params.json");
          param->Get<StochasticParams>()->LoadParams("params.json");
        }
      );
    auto* rm = simulation.GetResourceManager();

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
    std::uniform_real_distribution<real_t> offset(-physics->initial_cloud_radius, physics->initial_cloud_radius);

    for (int i = 0; i < physics->num_particles; ++i) 
    {
      Real3 position{offset(rng), offset(rng), offset(rng)};
      auto* particle = new ChargedAntiparticle(position, physics->species);

      if (interaction->coulomb_enabled) { particle->AddBehavior(new CoulombInteraction()); }
      particle->AddBehavior(new LorentzForce(&field));
      if (stochastic->annihilation_enabled) { particle->AddBehavior(new AnnihilationInteraction()); }

      rm->AddAgent(particle);
    }

    simulation.GetScheduler()->Simulate(physics->num_steps);

    rm->ForEachAgent([](Agent* agent) 
    {
      auto* particle = dynamic_cast<ChargedAntiparticle*>(agent);
      if (particle == nullptr) { return; }
      const auto& pos = particle->GetPosition();
      std::cout << "final position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
    });

    std::cout << "Simulation completed successfully!" << std::endl;
    return 0;
  }
}
#endif  // ANTIMATTER_KERNEL_H_