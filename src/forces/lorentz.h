#ifndef FORCES_LORENTZ_H_
#define FORCES_LORENTZ_H_

#include "core/behavior/behavior.h"
#include "core/container/math_array.h"
#include "core/param/param.h"
#include "core/util/math.h"
#include "core/simulation.h"
#include "params/physics_params.h"

namespace bdm
{
  class PenningTrapField;

  class LorentzForce : public Behavior
  {
    BDM_BEHAVIOR_HEADER(LorentzForce, Behavior, 1);

  public:
    LorentzForce() = default;
    explicit LorentzForce(const PenningTrapField* field)
        : field_(field),
          dt_(Simulation::GetActive()->GetParam()->Get<PhysicsParams>()->dt) {}
    ~LorentzForce() override = default;

    void Run(Agent* agent) override;

    const PenningTrapField* GetField() const { return field_; }
    void SetField(const PenningTrapField* field) { field_ = field; }

  private:
    const PenningTrapField* field_ = nullptr;
    real_t dt_ = 1e-9;
    bool velocity_staggered_ = false;  // true once v has been moved back to t = -dt/2
  };
}  // namespace bdm

#endif  // FORCES_LORENTZ_H_