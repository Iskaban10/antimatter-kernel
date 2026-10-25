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
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <chrono>

#include "biodynamo.h"

#include "agents/charged_antiparticle.h"
#include "analysis/energy_analysis.h"
#include "core/param/param.h"
#include "fields/penning_trap.h"
#include "forces/lorentz.h"
#include "interaction/annihilation.h"
#include "interaction/barnes-hut.h"
#include "interaction/coulomb_interaction.h"
#include "interaction/fast-multipole.h"
#include "params/interaction_params.h"
#include "params/physics_params.h"
#include "params/stochastic_params.h"

namespace bdm 
{
  inline int Simulate(int argc, const char** argv) 
  {
    auto start_time = std::chrono::high_resolution_clock::now();
    // Parameters that are to be loaded for simulation
    Param::RegisterParamGroup(new PhysicsParams());
    Param::RegisterParamGroup(new InteractionParams());
    Param::RegisterParamGroup(new StochasticParams());

    Simulation simulation(
        argc, argv, [](Param* param) {
          param->Get<PhysicsParams>()->LoadParams("params.json");
          param->Get<InteractionParams>()->LoadParams("params.json");
          param->Get<StochasticParams>()->LoadParams("params.json");
        });
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

    const CoulombMethod coulomb_method = interaction->GetCoulombMethod();
    const bool use_barnes_hut = interaction->coulomb_enabled && coulomb_method == CoulombMethod::kBarnesHut;
    const bool use_fast_multipole = interaction->coulomb_enabled && coulomb_method == CoulombMethod::kFastMultipole;

    const real_t energy_cutoff_radius = !interaction->coulomb_enabled ? static_cast<real_t>(0.0) : (coulomb_method == CoulombMethod::kStandard ? interaction->cutoff_radius : std::numeric_limits<real_t>::infinity());

    PenningTrapField field = PenningTrapField::MakeFromGeometry( physics->b0, physics->trap_voltage, physics->z0, physics->rho0);

    std::default_random_engine rng(physics->seed);
    std::uniform_real_distribution<real_t> offset(-physics->initial_cloud_radius, physics->initial_cloud_radius);

    for (const auto& [species, count] : physics->populations) 
    {
      for (int i = 0; i < count; i++) 
      {
        Real3 position{offset(rng), offset(rng), offset(rng)};
        auto* particle = new ChargedAntiparticle(position, species);

        if (interaction->coulomb_enabled && coulomb_method == CoulombMethod::kStandard)
          particle->AddBehavior(new CoulombInteraction());

        particle->AddBehavior(new LorentzForce(&field));

        if (stochastic->annihilation_enabled)   particle->AddBehavior(new AnnihilationInteraction());

        rm->AddAgent(particle);
      }
    }

    const int energy_log_interval = physics->num_steps >= 100 ? physics->num_steps / 100 : 1; // aprrox 1 % of number of steps

    const std::string energy_csv_path = simulation.GetOutputDir() + "/energy_analysis.csv";
    std::ofstream energy_csv(energy_csv_path);
    if (!energy_csv.is_open()) 
    {
      std::cerr << "WARNING: could not open '" << energy_csv_path << "' for writing -- energy CSV will not be saved this run." << std::endl;
    } 
    else 
    {
      energy_csv << "step,kinetic_energy_J,trap_potential_energy_J,coulomb_potential_energy_J,total_energy_J,num_active_particles,fractional_drift\n";
    }

    EnergyBreakdown initial_energy = ComputeTotalEnergy(rm, field, energy_cutoff_radius);
    std::cout << "\nstep 0  KE=" << initial_energy.kinetic
              << " PE_trap=" << initial_energy.trap_potential
              << " PE_coulomb=" << initial_energy.coulomb_potential
              << " TOTAL=" << initial_energy.Total() << " J" << std::endl;
    if (energy_csv.is_open()) 
    {
      energy_csv << 0 << "," << initial_energy.kinetic << ","
                << initial_energy.trap_potential << ","
                << initial_energy.coulomb_potential << ","
                << initial_energy.Total() << ","
                << initial_energy.num_active_particles << "," << 0.0 << "\n";
    }

    for (int step = 1; step <= physics->num_steps; ++step) 
    {
      if (use_barnes_hut)   ComputeBarnesHutCoulombForces(rm, interaction->barnes_hut_theta);
      else if (use_fast_multipole)  ComputeFastMultipoleCoulombForces(rm, interaction->fast_multipole_order, interaction->fast_multipole_theta);

      simulation.GetScheduler()->Simulate(1);

      if (step % energy_log_interval == 0 || step == physics->num_steps) 
      {
        EnergyBreakdown e = ComputeTotalEnergy(rm, field, energy_cutoff_radius);
        const real_t fractional_drift = initial_energy.Total() != 0.0
                ? (e.Total() - initial_energy.Total()) / std::abs(initial_energy.Total())
                : 0.0;
        std::cout << "step " << step << "  KE=" << e.kinetic
                  << " PE_trap=" << e.trap_potential
                  << " PE_coulomb=" << e.coulomb_potential
                  << " TOTAL=" << e.Total() << " J"
                  << "  (fractional drift from step 0: " << fractional_drift << ")" << std::endl;
        
        // write energy log results to the csv file
        if (energy_csv.is_open()) 
        {
          energy_csv << step << "," << e.kinetic << "," << e.trap_potential
                    << "," << e.coulomb_potential << "," << e.Total() << ","
                    << e.num_active_particles << "," << fractional_drift
                    << "\n";
        }
      }
    }

    if (energy_csv.is_open()) 
    {
      energy_csv.close();
      std::cout << "Energy analysis written to: " << energy_csv_path << std::endl;
    }

    rm->ForEachAgent([](Agent* agent) 
    {
      auto* particle = FastCast(agent);
      if (particle == nullptr)  return;
      const auto& pos = particle->GetPosition();
      std::cout << "final position: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
    });

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> simulation_duration = end_time - start_time;
    std::cout << "Simulation completed successfully!" << std::endl;
    std::cout << "Time taken "<< simulation_duration.count() << "s" << std::endl;
    return 0;
  }
}

#endif