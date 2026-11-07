#include "fields/penning_trap.h"

#include <cassert>
#include <cmath>
#include <random>

namespace bdm 
{
  PenningTrapField::PenningTrapField(real_t b0, real_t trap_voltage, real_t characteristic_dimension)
      : b0_(b0), trap_voltage_(trap_voltage), d_(characteristic_dimension) 
  {
    if (d_ <= 0.0) throw std::invalid_argument("characteristic dimension must be > 0");
  }

  PenningTrapField PenningTrapField::MakeFromGeometry(real_t b0, real_t trap_voltage, real_t z0, real_t rho0) 
  {
    const real_t d_squared = (z0 * z0) / 2.0 + (rho0 * rho0) / 4.0;
    const real_t d = std::sqrt(d_squared);
    return PenningTrapField(b0, trap_voltage, d);
  }

  Real3 PenningTrapField::GetElectricField(const Real3& position) const 
  {
    const real_t x = position[0];
    const real_t y = position[1];
    const real_t z = position[2];

    const real_t d_squared = d_ * d_;
    const real_t radial_coeff = trap_voltage_ / (2.0 * d_squared);
    const real_t axial_coeff = trap_voltage_ / d_squared;

    return Real3{radial_coeff * x, radial_coeff * y, -axial_coeff * z};
  }

  real_t PenningTrapField::GetPotential(const Real3& position) const 
  {
    const real_t x = position[0];
    const real_t y = position[1];
    const real_t z = position[2];
    const real_t rho_squared = x * x + y * y;
    const real_t d_squared = d_ * d_;
    return (trap_voltage_ / (4.0 * d_squared)) * (2.0 * z * z - rho_squared);
  }

  Real3 PenningTrapField::GetMagneticField(const Real3& /*position*/) const 
  {
    return Real3{0.0, 0.0, b0_};  // uniform axial field: independent of position for an ideal trap.
  }

  // set velocities for all the particles inside the trap
  std::vector<Real3> PenningTrapField::SetParticleVelocities(const std::map<ParticleSpecies, size_t>& populations, const std::map<ParticleSpecies, real_t>& particle_temperatures, const int seed)
  {
    size_t total = 0;
    for (const auto& entry : populations)   total += static_cast<size_t>(entry.second);
    if (total <= 0)   throw std::invalid_argument("PenningTrapField: Total particles in PenningTrap should be > 0");
  
    std::vector<Real3> velocities;
    velocities.reserve(total);
  
    for (const auto& [species, count] : populations) 
    {
      auto temp_it = particle_temperatures.find(species);
      double initial_temp = (temp_it != particle_temperatures.end()) ? temp_it->second : 0.0;

      if (initial_temp == 0.0) 
      {
        // since std::normal_distribution requires sigma > 0, so handle T = 0 as in a cold species is simply at rest.
        velocities.insert(velocities.end(), static_cast<size_t>(count), Real3{0.0, 0.0, 0.0});
        continue;
      }

      const double mass = static_cast<double>(GetDefaultSpeciesProperties(species).mass);
      // Thermal speed scale per Cartesian component, in m/s, depends on the species' mass.
      const double sigma = std::sqrt(static_cast<double>(constants::kBoltzmannConstant) * initial_temp / mass);

      std::random_device rd;
      std::mt19937 rng(seed);
      std::normal_distribution<double> gauss(0.0, sigma);
      std::vector<double> raw(3 * count);
      for (size_t i = 0; i < 3 * count; ++i)  raw[i] = gauss(rng);

      if (count >= 2) 
      {
        double mean[3] = {0.0, 0.0, 0.0};
        for (size_t i = 0; i < count; ++i) 
        {
          for (int d = 0; d < 3; ++d)   
          {
            mean[d] += raw[3 * i + d];
          }
        }
        for (int d = 0; d < 3; ++d)   
        {
          mean[d] /= static_cast<double>(count);
        }
        // Rescaling by sqrt(count/(count-1)) to restore the expected mean kinetic energy per particle which was lowered due to the constraint of mean drift = 0
        const double rescale = std::sqrt(static_cast<double>(count) / static_cast<double>(count - 1));
        for (size_t i = 0; i < count; ++i) 
        {
          for (int d = 0; d < 3; ++d)     
          {
            raw[3 * i + d] = (raw[3 * i + d] - mean[d]) * rescale;
          }
        }
      }
  
      for (size_t i = 0; i < count; ++i) 
      {
        velocities.push_back(Real3{static_cast<real_t>(raw[3 * i]),
                                  static_cast<real_t>(raw[3 * i + 1]),
                                  static_cast<real_t>(raw[3 * i + 2])});
      }
    }
    return velocities;
  }

  // set positions for all the particles inside the trap
  std::vector<Real3> PenningTrapField::SetParticlePositions(const std::map<ParticleSpecies, size_t>& populations, const Real3& radii, const int seed)
  {
    std::vector<Real3> positions;
  
    // Calculate total number of particles to pre-allocate vector memory
    size_t total_particles = 0;
    for (const auto& [species, count] : populations)    
      total_particles += count;
    positions.reserve(total_particles);

    // Initialize the random number engine
    std::default_random_engine rng(seed);
    std::uniform_real_distribution<real_t> offset_x(-radii[0], radii[0]);
    std::uniform_real_distribution<real_t> offset_y(-radii[1], radii[1]);
    std::uniform_real_distribution<real_t> offset_z(-radii[2], radii[2]);
  
    // Generate coordinates for each particle category block
  for (const auto& [species, count] : populations) 
  {
    for (size_t i = 0; i < count; ++i) 
    {
      Real3 point;
      bool inside_ellipsoid = false;

      // Rejection sampling loop
      while (!inside_ellipsoid) 
      {
        point[0] = offset_x(rng);
        point[1] = offset_y(rng);
        point[2] = offset_z(rng);

        // Ellipsoid normalization math layout: (x/a)² + (y/b)² + (z/c)² <= 1
        real_t norm_x = point[0] / radii[0];
        real_t norm_y = point[1] / radii[1];
        real_t norm_z = point[2] / radii[2];

        if ((norm_x * norm_x) + (norm_y * norm_y) + (norm_z * norm_z) <= 1.0) 
          inside_ellipsoid = true; // Coordinate verified
      }
      positions.push_back(point);
    }
  }
    return positions;
  }
}