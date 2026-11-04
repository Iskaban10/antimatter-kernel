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
        case ParticleSpecies::kAntihelium: return "antihelium";
        default: return "unknown";
      }
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
        populations[ParseSpecies(name)] = count.get<size_t>();
      }
    }

    if (j.contains("particle_temperatures")) 
    {
      particle_temperatures.clear();
      for (const auto& [name, count] : j.at("particle_temperatures").items()) 
      {
        particle_temperatures[ParseSpecies(name)] = count.get<size_t>();
      }
    }

    num_steps = j.value("num_steps", num_steps);
    if (j.contains("initial_cloud_radii")) 
    {
      const auto& radii_json = j.at("initial_cloud_radii");
      
      if (radii_json.contains("radius_x")) {
        initial_cloud_radii[0] = radii_json.at("radius_x").get<real_t>();
      }
      if (radii_json.contains("radius_y")) {
        initial_cloud_radii[1] = radii_json.at("radius_y").get<real_t>();
      }
      if (radii_json.contains("radius_z")) {
        initial_cloud_radii[2] = radii_json.at("radius_z").get<real_t>();
      }
    }
  }

  void PhysicsParams::PrintParams() const 
  {
    std::cout << "[physics] \nseed = " << seed
              << "\nb0 = " << b0 << " T"
              << "\ntrap_voltage = " << trap_voltage << " V"
              << "\nz0 = " << z0 << " m"
              << "\nrho0 = " << rho0 << " m"
              << "\ndt = " << dt << " s"
              << "\nnum_steps = " << num_steps << std::endl;

    std::cout << "\nPARTICLE POPULATIONS:";
    for (const auto& [species, count] : populations) 
    {
      std::cout << '\n' << SpeciesToString(species) << " x " << count;
    }
    std::cout<<'\n';
    std::cout << "\nPARTICLE TEMPERATURES:";
    for (const auto& [species, temp] : particle_temperatures)
    {
      std::cout << '\n' << SpeciesToString(species) << " - " << temp << " Kelvin";
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
    for (const auto& [species, temperature] : particle_temperatures) 
    {
      if (temperature < 0)
        throw std::invalid_argument("PhysicsParams: initial temperature of particles count must be >= 0");
    }
    if (populations.empty())
      throw std::invalid_argument("PhysicsParams: populations must not be empty");
    for (const auto& [species, count] : populations) 
    {
      if (count < 0)
        throw std::invalid_argument("PhysicsParams: every population count must be >= 0");
    }
    // initial cloud radius checks
    if (initial_cloud_radii[0] <= 0.0 || initial_cloud_radii[1] <= 0.0 || initial_cloud_radii[2] <= 0.0) 
    {
      throw std::invalid_argument("PhysicsParams: All components of initial_cloud_radii must be > 0.");
    }
    if (initial_cloud_radii[0] >= rho0 || initial_cloud_radii[1] >= rho0) 
    {
      throw std::invalid_argument("PhysicsParams: Radial radii (X and Y) must be well inside trap radius (rho0).");
    }
    if (initial_cloud_radii[2] >= z0) 
    {
      throw std::invalid_argument("PhysicsParams: Axial radius (Z) must be well inside trap length (z0).");
    }
  }
}