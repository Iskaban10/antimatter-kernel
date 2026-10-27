#include "params/interaction_params.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "json.hpp"

namespace bdm 
{

  const ParamGroupUid InteractionParams::kUid = ParamGroupUidGenerator::Get()->NewUid();

  ComputationMethod InteractionParams::GetCoulombMethod() const 
  {
    if (coulomb_method == "standard") return ComputationMethod::kStandard;
    if (coulomb_method == "barnes-hut") return ComputationMethod::kBarnesHut;
    if (coulomb_method == "fast-multipole") return ComputationMethod::kFastMultipole;
    throw std::invalid_argument( "InteractionParams: unrecognized coulomb_method '" + coulomb_method + "'. Expected one of: standard, barnes-hut, fast-multipole.");
  }

  ComputationMethod InteractionParams::GetEnergyMethod() const 
  {
    if (energy_method == "standard") return ComputationMethod::kStandard;
    if (energy_method == "barnes-hut") return ComputationMethod::kBarnesHut;
    if (energy_method == "fast-multipole") return ComputationMethod::kFastMultipole;
    throw std::invalid_argument( "InteractionParams: unrecognized energy_method '" + energy_method + "'. Expected one of: standard, barnes-hut, fast-multipole.");
  }

  void InteractionParams::LoadParams(const std::string& filename) 
  {
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
      std::cout << "InteractionParams::LoadParams: '" << filename << "' not found -- using default values." << std::endl;
      return;
    }
    nlohmann::json j;
    file >> j;
    coulomb_enabled = j.value("coulomb_enabled", coulomb_enabled);
    coulomb_method = j.value("coulomb_computation_method", coulomb_method);
    energy_method = j.value("energy_computation_method", energy_method);
    cutoff_radius = j.value("cutoff_radius", cutoff_radius);
    barnes_hut_theta = j.value("barnes_hut_theta", barnes_hut_theta);
    fast_multipole_order = j.value("fast_multipole_order", fast_multipole_order);
    fast_multipole_theta = j.value("fast_multipole_theta", fast_multipole_theta);
  }

  void InteractionParams::PrintParams() const 
  {
    std::cout << "[interaction] \ncoulomb_enabled=" << (coulomb_enabled ? "true" : "false") << "\ncoulomb_computation_method=" << coulomb_method;
    if (coulomb_enabled)
    {
      if (coulomb_method == "standard") 
      {
        std::cout << "\ncutoff_radius=" << cutoff_radius << " m";
      } 
      else if (coulomb_method == "barnes-hut") 
      {
        std::cout << "\nbarnes_hut_theta=" << barnes_hut_theta;
      } 
      else if (coulomb_method == "fast-multipole") 
      {
        std::cout << "\nfast_multipole_order=" << fast_multipole_order << "\nfast_multipole_theta=" << fast_multipole_theta;
      }
    }
    std::cout << "\nenergy_computation_method=" << energy_method;
    std::cout << std::endl;
  }

  void InteractionParams::Validate() const 
  {
    const ComputationMethod coulomb_method = GetCoulombMethod(); 
    const ComputationMethod energy_method = GetCoulombMethod(); 

    if (energy_method == ComputationMethod::kStandard && cutoff_radius <= 0.0) 
    {
      throw std::invalid_argument( "InteractionParams: cutoff_radius must be > 0 when energy_method is 'standard' and coulomb_enabled is true.");
    }
    if (energy_method == ComputationMethod::kBarnesHut && (barnes_hut_theta < 0.0 || barnes_hut_theta > 1.0)) 
    {
      throw std::invalid_argument( "InteractionParams: barnes_hut_theta must be in [0, 1] (values above 1 could make a particle interact with its own tree node).");
    }
    if (energy_method == ComputationMethod::kFastMultipole) 
    {
      if (fast_multipole_order < 1 || fast_multipole_order > 8) 
      {
        throw std::invalid_argument("InteractionParams: fast_multipole_order must be an integer in [1, 8].");
      }
      if (fast_multipole_theta < 0.0 || fast_multipole_theta >= 1.0) 
      {
        throw std::invalid_argument("InteractionParams: fast_multipole_theta must be in [0, 1) (at 1 or above the expansions no longer converge).");
      }
    }
    if (!coulomb_enabled)   return;
    if (coulomb_method == ComputationMethod::kStandard && cutoff_radius <= 0.0) 
    {
      throw std::invalid_argument( "InteractionParams: cutoff_radius must be > 0 when coulomb_method " "is 'standard' and coulomb_enabled is true.");
    }
    if (coulomb_method == ComputationMethod::kBarnesHut && (barnes_hut_theta < 0.0 || barnes_hut_theta > 1.0)) 
    {
      throw std::invalid_argument( "InteractionParams: barnes_hut_theta must be in [0, 1] (values " "above 1 could make a particle interact with its own tree node).");
    }
    if (coulomb_method == ComputationMethod::kFastMultipole) 
    {
      if (fast_multipole_order < 1 || fast_multipole_order > 8) 
      {
        throw std::invalid_argument("InteractionParams: fast_multipole_order must be an integer in [1, 8].");
      }
      if (fast_multipole_theta < 0.0 || fast_multipole_theta >= 1.0) 
      {
        throw std::invalid_argument("InteractionParams: fast_multipole_theta must be in [0, 1) (at 1 or above " "the expansions no longer converge).");
      }
    }
  }
}