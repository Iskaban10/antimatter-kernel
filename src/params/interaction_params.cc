#include "params/interaction_params.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "json.hpp"

namespace bdm 
{
  const ParamGroupUid InteractionParams::kUid = ParamGroupUidGenerator::Get()->NewUid();
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
    cutoff_radius = j.value("cutoff_radius", cutoff_radius);
  }

  void InteractionParams::PrintParams() const 
  {
    std::cout << "[interaction] coulomb_enabled=" << (coulomb_enabled ? "true" : "false") << " cutoff_radius=" << cutoff_radius << " m (PLACEHOLDER -- not physically derived)" << std::endl;
  }

  void InteractionParams::Validate() const 
  {
    if (coulomb_enabled && cutoff_radius <= 0.0) 
    {
      throw std::invalid_argument(
          "InteractionParams: cutoff_radius must be > 0 when coulomb_enabled "
          "is true.");
    }
  }
}