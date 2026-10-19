#include "params/stochastic_params.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "json.hpp"

namespace bdm 
{
  const ParamGroupUid StochasticParams::kUid = ParamGroupUidGenerator::Get()->NewUid();

  void StochasticParams::LoadParams(const std::string& filename) 
  {
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
      std::cout << "StochasticParams::LoadParams: '" << filename << "' not found -- using default values." << std::endl;
      return;
    }
    nlohmann::json j;
    file >> j;
    annihilation_enabled = j.value("annihilation_enabled", annihilation_enabled);
    contact_radius = j.value("contact_radius", contact_radius);
    annihilation_probability = j.value("annihilation_probability", annihilation_probability);
  }

  void StochasticParams::PrintParams() const 
  {
    std::cout << "\n[stochastic] \nannihilation_enabled=" << (annihilation_enabled ? "true" : "false")
              << " \ncontact_radius=" << contact_radius << " m"
              << " \nannihilation_probability=" << annihilation_probability << std::endl;
  }

  void StochasticParams::Validate() const 
  {
    if (!annihilation_enabled) return;
    if (contact_radius <= 0.0) 
    {
      throw std::invalid_argument( "StochasticParams: contact_radius must be > 0 when " "annihilation_enabled is true.");
    }
    if (annihilation_probability < 0.0 || annihilation_probability > 1.0) 
    {
      throw std::invalid_argument("StochasticParams: annihilation_probability must be in [0, 1].");
    }
  }
}