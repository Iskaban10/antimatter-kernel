#include "params/interaction_params.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "json.hpp"

namespace bdm 
{

  const ParamGroupUid InteractionParams::kUid = ParamGroupUidGenerator::Get()->NewUid();

  CoulombMethod InteractionParams::GetCoulombMethod() const 
  {
    if (coulomb_method == "standard") return CoulombMethod::kStandard;
    if (coulomb_method == "barnes-hut") return CoulombMethod::kBarnesHut;
    if (coulomb_method == "fast-multipole") return CoulombMethod::kFastMultipole;
    throw std::invalid_argument( "InteractionParams: unrecognized coulomb_method '" + coulomb_method + "'. Expected one of: standard, barnes-hut, fast-multipole.");
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
    coulomb_method = j.value("coulomb_method", coulomb_method);
    cutoff_radius = j.value("cutoff_radius", cutoff_radius);
    barnes_hut_theta = j.value("barnes_hut_theta", barnes_hut_theta);
    fast_multipole_order = j.value("fast_multipole_order", fast_multipole_order);
    fast_multipole_theta = j.value("fast_multipole_theta", fast_multipole_theta);
  }

  void InteractionParams::PrintParams() const 
  {
    std::cout << "\n[interaction] coulomb_enabled=" << (coulomb_enabled ? "true" : "false") << "\ncoulomb_method=" << coulomb_method;
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
    std::cout << std::endl;
  }

  void InteractionParams::Validate() const 
  {
    const CoulombMethod method = GetCoulombMethod();  // throws on a bad name

    if (!coulomb_enabled)   return;
    if (method == CoulombMethod::kStandard && cutoff_radius <= 0.0) 
    {
      throw std::invalid_argument( "InteractionParams: cutoff_radius must be > 0 when coulomb_method " "is 'standard' and coulomb_enabled is true.");
    }
    if (method == CoulombMethod::kBarnesHut && (barnes_hut_theta < 0.0 || barnes_hut_theta > 1.0)) 
    {
      throw std::invalid_argument( "InteractionParams: barnes_hut_theta must be in [0, 1] (values " "above 1 could make a particle interact with its own tree node).");
    }
    if (method == CoulombMethod::kFastMultipole) 
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