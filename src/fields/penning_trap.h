#ifndef FIELDS_PENNING_TRAP_H_
#define FIELDS_PENNING_TRAP_H_

#include "core/container/math_array.h"

namespace bdm 
{
  class PenningTrapField 
  {
    public:
      PenningTrapField() = default;

      // constructs a field directly from B0, trap voltage U0, and the precomputed characteristic trap-dimension d.
      PenningTrapField(real_t b0, real_t trap_voltage, real_t characteristic_dimension);

      // constructs a field from physical electrode geometry and computes the characteristic dimension d internally.
      static PenningTrapField MakeFromGeometry(real_t b0, real_t trap_voltage, real_t z0, real_t rho0);

      // Returns the electric field E(x, y, z) at `position`, in V/m.
      Real3 GetElectricField(const Real3& position) const;
      Real3 GetMagneticField(const Real3& position) const;

      real_t GetB0() const { return b0_; }
      void SetB0(real_t b0) { b0_ = b0; }

      real_t GetTrapVoltage() const { return trap_voltage_; }
      void SetTrapVoltage(real_t trap_voltage) { trap_voltage_ = trap_voltage; }

      real_t GetCharacteristicDimension() const { return d_; }
      void SetCharacteristicDimension(real_t d) { d_ = d; }

    private:
      // axial magnetic field strength, Tesla.
      real_t b0_ = 0.0;

      // trap (ring-to-endcap) voltage, Volts.
      real_t trap_voltage_ = 0.0;

      // characteristic trap dimension, meters. Must be > 0 for GetElectricField to return finite values
      real_t d_ = 1.0;
  };
}
#endif