#ifndef FIELDS_PENNING_TRAP_H_
#define FIELDS_PENNING_TRAP_H_

#include <vector>
#include <map>
#include "agents/particle_species.h"
#include "core/container/math_array.h"

namespace bdm 
{

  // Describes the static electromagnetic field configuration of an ideal (hyperbolic-electrode) Penning trap:
  // a uniform axial magnetic field superimposed on a quadrupole electrostatic confining potential.
  // This class is deliberately field-geometry-only answering only the E and B in space
  //
  // Field model:
  //   B(x, y, z) = B0 * z_hat  (uniform)
  //
  //   phi(rho, z) = (U0 / (4 d^2)) * (2 z^2 - rho^2),  rho^2 = x^2 + y^2
  //   E = -grad(phi)
  //   Ex = (U0 / (2 d^2)) * x
  //   Ey = (U0 / (2 d^2)) * y
  //   Ez = -(U0 / d^2) * z
  //
  // where d is the trap's characteristic dimension, d^2 = z0^2/2 + rho0^2/4, with z0 the endcap half-distance and rho0 the ring-electrode radius.

  class PenningTrapField 
  {
    public:
      PenningTrapField() = default;

      // Constructs a field directly from B0, trap voltage U0, and the precomputed characteristic dimension d.
      PenningTrapField(real_t b0, real_t trap_voltage, real_t characteristic_dimension);

      // Constructs a field from physical electrode geometry (endcap half-distance z0 and ring-electrode radius rho0),
      // and computes the characteristic dimension d internally.
      static PenningTrapField MakeFromGeometry(real_t b0, real_t trap_voltage, real_t z0, real_t rho0);

      // Returns the electric field E(x, y, z) at `position`, in V/m.
      Real3 GetElectricField(const Real3& position) const;

      // Returns the electric potential phi(x, y, z) at `position`, in Volts (phi = (U0/4d^2)(2z^2 - rho^2), matching GetElectricField = -grad(phi)).
      real_t GetPotential(const Real3& position) const;

      // Returns the magnetic field B(x, y, z) at `position`, in Tesla.
      Real3 GetMagneticField(const Real3& position) const; // position unrequired in case of uniform field

      real_t GetB0() const { return b0_; }
      void SetB0(real_t b0) { b0_ = b0; }

      real_t GetTrapVoltage() const { return trap_voltage_; }
      void SetTrapVoltage(real_t trap_voltage) { trap_voltage_ = trap_voltage; }

      real_t GetCharacteristicDimension() const { return d_; }
      void SetCharacteristicDimension(real_t d) { d_ = d; }

      static std::vector<Real3> SetParticleVelocities(const std::map<ParticleSpecies, size_t>& populations, real_t initial_temp);

    private:
      // Axial magnetic field strength, Tesla.
      real_t b0_ = 0.0;

      // Trap (ring-to-endcap) voltage, Volts.
      real_t trap_voltage_ = 0.0;

      // Characteristic trap dimension, meters. Must be > 0 for
      // GetElectricField to return finite values; see MakeFromGeometry.
      real_t d_ = 1.0;
  };
}
#endif