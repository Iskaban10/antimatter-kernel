#include "params/load_params.h"

#include <fstream>
#include <iostream>
#include "json.hpp"

namespace bdm 
{
  namespace 
  {
    // converts a species name string from params.json into a ParticleSpecies. 
    ParticleSpecies ParseSpecies(const std::string& name) 
    {
      if (name == "antiproton") return ParticleSpecies::kAntiproton;
      if (name == "positron") return ParticleSpecies::kPositron;
      if (name == "antihydrogen") return ParticleSpecies::kAntihydrogen;
      if (name == "electron") return ParticleSpecies::kElectron;
      if (name == "proton") return ParticleSpecies::kProton;
      throw std::invalid_argument( "load_params: unrecognized species name '" + name +
          "' in params.json. Expected one of: antiproton, positron, "
          "antihydrogen, electron, proton.");
    }

    // converts particle species to string for printing
    std::string SpeciesToString(ParticleSpecies species) 
    {
      switch (species) 
      {
        case ParticleSpecies::kAntiproton: return "antiproton";
        case ParticleSpecies::kPositron: return "positron";
        case ParticleSpecies::kAntihydrogen: return "antihydrogen";
        case ParticleSpecies::kElectron: return "electron";
        case ParticleSpecies::kProton: return "proton";
      }
      return "unknown";
    }
  }

  void LoadParams(const std::string& path, PhysicsParams* physics, InteractionParams* interaction, StochasticParams* stochastic) 
  {
    std::ifstream file(path);
    if (!file.is_open()) 
    {
      std::cout << "load_params: '" << path << "' not found - using default parameter values." << std::endl;
      return;
    }

    nlohmann::json j;
    file >> j;

    physics->seed = j.value("seed", physics->seed);
    physics->b0 = j.value("b0", physics->b0);
    physics->trap_voltage = j.value("trap_voltage", physics->trap_voltage);
    physics->z0 = j.value("z0", physics->z0);
    physics->rho0 = j.value("rho0", physics->rho0);
    physics->dt = j.value("dt", physics->dt);
    physics->num_steps = j.value("num_steps", physics->num_steps);
    physics->num_particles = j.value("num_particles", physics->num_particles);
    physics->initial_cloud_radius = j.value("initial_cloud_radius", physics->initial_cloud_radius);

    if (j.contains("species")) { physics->species = ParseSpecies(j.at("species").get<std::string>());}

    interaction->coulomb_enabled = j.value("coulomb_enabled", interaction->coulomb_enabled);
    interaction->cutoff_radius = j.value("cutoff_radius", interaction->cutoff_radius);

    stochastic->annihilation_enabled = j.value("annihilation_enabled", stochastic->annihilation_enabled);
    stochastic->contact_radius = j.value("contact_radius", stochastic->contact_radius);
    stochastic->annihilation_probability = j.value( "annihilation_probability", stochastic->annihilation_probability);
    std::cout << "load_params: loaded configuration from '" << path << "'." << std::endl;
  }

  void PrintParams(const PhysicsParams& physics, const InteractionParams& interaction, const StochasticParams& stochastic) 
  {
      std::cout << "------------- AntimatterKernel Parameters -------------" << std::endl;
      std::cout << "[physics]" << std::endl;
      std::cout << "  seed                 = " << physics.seed << std::endl;
      std::cout << "  species              = " << SpeciesToString(physics.species) << std::endl;
      std::cout << "  num_particles        = " << physics.num_particles << std::endl;
      std::cout << "  initial_cloud_radius = " << physics.initial_cloud_radius << " m" << std::endl;
      std::cout << "  b0                   = " << physics.b0 << " T" << std::endl;
      std::cout << "  trap_voltage         = " << physics.trap_voltage << " V" << std::endl;
      std::cout << "  z0                   = " << physics.z0 << " m" << std::endl;
      std::cout << "  rho0                 = " << physics.rho0 << " m" << std::endl;
      std::cout << "  dt                   = " << physics.dt << " s" << std::endl;
      std::cout << "  num_steps            = " << physics.num_steps << std::endl;

      std::cout << "[interaction]" << std::endl;
      std::cout << "  coulomb_enabled      = " << (interaction.coulomb_enabled ? "true" : "false") << std::endl;
      std::cout << "  cutoff_radius        = " << interaction.cutoff_radius << " m" << std::endl;

      std::cout << "[stochastic]" << std::endl;
      std::cout << "  annihilation_enabled = " << (stochastic.annihilation_enabled ? "true" : "false") << std::endl;
      std::cout << "  contact_radius       = " << stochastic.contact_radius << " m" << std::endl;
      std::cout << "  annihilation_probability = " << stochastic.annihilation_probability << std::endl;
      std::cout << "-------------------------------------------" << std::endl;
  }
}  // namespace bdm