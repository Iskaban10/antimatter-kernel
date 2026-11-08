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
#include "antimatter-kernel.h"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "biodynamo.h"

#include "agents/charged_antiparticle.h"
#include "analysis/energy_analysis.h"
#include "core/param/param.h"
#include "fields/penning_trap.h"
#include "forces/cooling.h"
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
  int Simulate(int argc, const char** argv) 
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
    
    std::cout<<std::endl;
    physics->PrintParams();
    std::cout<<std::endl;
    interaction->PrintParams();
    std::cout<<std::endl;
    stochastic->PrintParams();

    // Cyclotron cooling removes energy on purpose, so it must be OFF when
    // validating energy conservation of the integrator. Set to true for real physics runs.

    const ComputationMethod coulomb_method = interaction->GetCoulombMethod();
    const bool use_standard = interaction->coulomb_enabled && coulomb_method == ComputationMethod::kStandard;
    const bool use_barnes_hut = interaction->coulomb_enabled && coulomb_method == ComputationMethod::kBarnesHut;
    const bool use_fast_multipole = interaction->coulomb_enabled && coulomb_method == ComputationMethod::kFastMultipole;
    const real_t energy_cutoff_radius = interaction->cutoff_radius;

    PenningTrapField field = PenningTrapField::MakeFromGeometry( physics->b0, physics->trap_voltage, physics->z0, physics->rho0);
    const std::vector<Real3> particle_positions = PenningTrapField::SetParticlePositions( physics->populations, physics->initial_cloud_radii, physics->seed);  // array of random distributions of particles inside the initial cloud in the Penning trap
    const std::vector<Real3> particle_velocities = PenningTrapField::SetParticleVelocities( physics->populations, physics->particle_temperatures, physics->seed);  // array of gaussian distribution of particle velocities

    size_t particle_counter = 0; // counter to keep track on particles
    for (const auto& [species, count] : physics->populations) 
    {
      for (size_t i = 0; i < count; i++)
      {
        auto* particle = new ChargedAntiparticle(particle_positions[particle_counter], species, particle_velocities[particle_counter]);
        particle_counter++;
        particle->AddBehavior(new LorentzForce(&field));
        particle->AddBehavior(new CyclotronCooling(&field, /*tau_scale=*/physics->cooling_scale, /*bath_temp_K=*/physics->bath_temp));

        if (stochastic->annihilation_enabled)   particle->AddBehavior(new AnnihilationInteraction());

        rm->AddAgent(particle);
      }
    }

    // Energy Analysis
    const ComputationMethod energy_method = interaction->GetEnergyMethod();
    EnergyComputationSettings energy_settings = {energy_method, energy_cutoff_radius,
                                                    interaction->barnes_hut_theta, interaction->fast_multipole_theta,
                                                    interaction->fast_multipole_order};

    const size_t energy_log_interval = physics->num_steps >= 100 ? physics->num_steps / 100 : 1; // aprrox 1 % of number of steps

    // Time-centering needs the previous step's potentials, so this lives for the whole run.
    EnergyHistory energy_history;

    auto evaluate_energy = [&]() -> EnergyBreakdown
    {
      EnergyBreakdown e = ComputeTotalEnergy(rm, field, energy_settings, &energy_history);
      if (!interaction->coulomb_enabled)  e.energy_potential = 0.0;
      return e;
    };

    const std::string energy_csv_path = simulation.GetOutputDir() + "/energy_analysis_data.csv";
    std::ofstream energy_csv(energy_csv_path);
    if (!energy_csv.is_open())
    {
      std::cerr << "WARNING: could not open '" << energy_csv_path << "' for writing -- energy CSV will not be saved this run." << std::endl;
    } 
    else 
    {
      energy_csv << "step,kinetic_energy_J,trap_potential_energy_J,potential_energy_J,total_energy_J,num_active_particles,fractional_drift,absolute_drift_J\n";
    }

    EnergyBreakdown initial_energy = evaluate_energy();

    // Characteristic energy for the drift: the trap depth (every species carries charge +/- e).
    // Dividing by the initial total is unreliable because it can be tiny for cold, bound particles.
    const real_t elementary_charge = 1.602176634e-19;
    real_t energy_scale = elementary_charge * std::abs(physics->trap_voltage);
    if (energy_scale <= 0.0)  energy_scale = std::abs(initial_energy.Total());
    if (energy_scale <= 0.0)  energy_scale = 1.0;

    std::cout << "\nstep 0  KE=" << initial_energy.kinetic
              << " PE_trap=" << initial_energy.trap_potential
              << " PE_coulomb=" << initial_energy.energy_potential
              << " TOTAL=" << initial_energy.Total() << " J" << std::endl;
    if (energy_csv.is_open())
    {
      energy_csv << 0 << "," << initial_energy.kinetic << ","
                << initial_energy.trap_potential << ","
                << initial_energy.energy_potential << ","
                << initial_energy.Total() << ","
                << initial_energy.num_active_particles << "," << 0.0 << "," << 0.0 << "\n";
    }

    for (size_t step = 1; step <= physics->num_steps; ++step) 
    {
      if (use_barnes_hut)   
        ComputeBarnesHutCoulombForces(rm, interaction->barnes_hut_theta);
      else if (use_fast_multipole)  
        ComputeFastMultipoleCoulombForces(rm, interaction->fast_multipole_order, interaction->fast_multipole_theta);
      else if (use_standard)
        ComputeStandardCoulombForces(rm, interaction->cutoff_radius);

      simulation.GetScheduler()->Simulate(1);

      const bool log_step = (step % energy_log_interval == 0) || (step == physics->num_steps);

      // The time-centered energy pairs step n with step n-1, so the step just before
      // a log step must be evaluated too (even though it is not printed).
      const bool need_energy = log_step || ((step + 1) % energy_log_interval == 0) || (step + 1 == physics->num_steps);

      if (need_energy)
      {
        EnergyBreakdown e = evaluate_energy();

        if (log_step)
        {
          const real_t absolute_drift = e.Total() - initial_energy.Total();
          const real_t fractional_drift = absolute_drift / energy_scale;

          std::cout << "step " << step << "  KE=" << e.kinetic
                    << " PE_trap=" << e.trap_potential
                    << " PE_coulomb=" << e.energy_potential
                    << " TOTAL=" << e.Total() << " J"
                    << "  (drift: " << absolute_drift << " J, "
                    << fractional_drift << " of trap depth)" << std::endl;

          // write energy log results to the csv file
          if (energy_csv.is_open()) 
          {
            energy_csv << step << "," << e.kinetic << "," << e.trap_potential
                      << "," << e.energy_potential << "," << e.Total() << ","
                      << e.num_active_particles << "," << fractional_drift << ","
                      << absolute_drift << "\n";
          }
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
      auto* particle = static_cast<ChargedAntiparticle*>(agent);
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
}  // namespace bdm

int main(int argc, const char** argv) { return bdm::Simulate(argc, argv); }