#include "fields/penning_trap.h"

#include <cassert>
#include <cmath>

namespace bdm 
{
  PenningTrapField::PenningTrapField(real_t b0, real_t trap_voltage, real_t characteristic_dimension) : b0_(b0), trap_voltage_(trap_voltage), d_(characteristic_dimension) 
  {
    assert(d_ > 0.0 && "Penning trap characteristic dimension must be > 0");
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

  Real3 PenningTrapField::GetMagneticField(const Real3& /*position*/) const 
  {
    return Real3{0.0, 0.0, b0_}; // uniform axial field independent of trap
  }
}  // namespace bdm