#include "params/physics_params.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "json.hpp"

namespace bdm 
{
  const ParamGroupUid PhysicsParams::kUid = ParamGroupUidGenerator::Get()->NewUid();
  namespace 
  {
    ParticleSpecies ParseSpecies(const std::string& name) 
    {
      if (name == "antiproton") return ParticleSpecies::kAntiproton;
      if (name == "positron") return ParticleSpecies::kPositron;
      if (name == "antihydrogen") return ParticleSpecies::kAntihydrogen;
      if (name == "electron") return ParticleSpecies::kElectron;
      if (name == "proton") return ParticleSpecies::kProton;
      throw std::invalid_argument("PhysicsParams::LoadParams: unrecognized species name '" + name + "'");
    }

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

  void PhysicsParams::LoadParams(const std::string& filename) 
  {
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
      std::cout << "PhysicsParams::LoadParams: '" << filename << "' not found -- using default values." << std::endl;
      return;
    }
    nlohmann::json j;
    file >> j;

    seed = j.value("seed", seed);
    b0 = j.value("b0", b0);
    trap_voltage = j.value("trap_voltage", trap_voltage);
    z0 = j.value("z0", z0);
    rho0 = j.value("rho0", rho0);
    dt = j.value("dt", dt);

    if (j.contains("particle_populations")) 
    {
      populations.clear();
      for (const auto& [name, count] : j.at("particle_populations").items()) 
      {
        populations[ParseSpecies(name)] = count.get<int>();
      }
    }

    num_steps = j.value("num_steps", num_steps);
    initial_cloud_radius = j.value("initial_cloud_radius", initial_cloud_radius);
  }

  void PhysicsParams::PrintParams() const 
  {
    std::cout << "\n[physics] \nseed = " << seed
              << "\ninitial_cloud_radius = " << initial_cloud_radius << " m"
              << "\nb0 = " << b0 << " T"
              << "\ntrap_voltage = " << trap_voltage << " V"
              << "\nz0 = " << z0 << " m"
              << "\nrho0 = " << rho0 << " m"
              << "\ndt = " << dt << " s"
              << "\nnum_steps = " << num_steps << std::endl;

    std::cout << "Particle populations:";
    for (const auto& [species, count] : populations) 
    {
      std::cout << "\t(" << SpeciesToString(species) << " x " << count << ")";
    }
    std::cout<<'\n';
  }

  void PhysicsParams::Validate() const 
  {
    if (b0 <= 0.0) throw std::invalid_argument("PhysicsParams: b0 must be > 0");
    if (z0 <= 0.0 || rho0 <= 0.0)
      throw std::invalid_argument("PhysicsParams: z0 and rho0 must be > 0");
    if (dt <= 0.0) throw std::invalid_argument("PhysicsParams: dt must be > 0");
    if (num_steps <= 0)
      throw std::invalid_argument("PhysicsParams: num_steps must be > 0");
    if (populations.empty())
      throw std::invalid_argument("PhysicsParams: populations must not be empty");
    for (const auto& [species, count] : populations) 
    {
      if (count < 0)
        throw std::invalid_argument("PhysicsParams: every population count must be > 0");
    }
    if (initial_cloud_radius <= 0.0 || initial_cloud_radius >= z0 || initial_cloud_radius >= rho0) 
    {
      throw std::invalid_argument( "PhysicsParams: initial_cloud_radius must be > 0 and well inside " "(z0, rho0)");
    }
  }
}