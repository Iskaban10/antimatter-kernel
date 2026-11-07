#ifndef FORCES_COOLING_H_
#define FORCES_COOLING_H_

#include "core/behavior/behavior.h"
#include "core/container/math_array.h"

namespace bdm 
{
  class PenningTrapField;

  // Cyclotron (synchrotron) radiation cooling.
  // Damps the velocity component perpendicular to B; the parallel component is untouched.
  //   1/tau_E = q^4 B^2 / (3 pi eps0 m^3 c^3)   (perpendicular energy e-folding time)
  //   v_perp decays as exp(-dt / (2 tau_E))
  // With bath_temp_K > 0, Langevin noise makes v_perp relax to a thermal floor instead of zero.
  class CyclotronCooling : public Behavior 
  {
    BDM_BEHAVIOR_HEADER(CyclotronCooling, Behavior, 1);

  public:
    CyclotronCooling() = default;
    // tau_scale   : multiplies tau_E (< 1 accelerates cooling to simulable timescales)
    // bath_temp_K : bath temperature in kelvin (0 = pure damping)
    explicit CyclotronCooling(const PenningTrapField* field, real_t tau_scale = 1.0, real_t bath_temp_K = 0.0)
        : field_(field), tau_scale_(tau_scale), bath_temp_K_(bath_temp_K) {}
    ~CyclotronCooling() override = default;

    void Run(Agent* agent) override;

    // Perpendicular-energy e-folding time (SI). Returns +inf if the rate is 0.
    static real_t CoolingTime(real_t charge, real_t mass, real_t b_mag);

    const PenningTrapField* GetField() const { return field_; }
    void SetField(const PenningTrapField* field) { field_ = field; }

  private:
    const PenningTrapField* field_ = nullptr;
    real_t tau_scale_ = 1.0;
    real_t bath_temp_K_ = 0.0;
  };
}

#endif