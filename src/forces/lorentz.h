#ifndef FORCES_LORENTZ_H_
#define FORCES_LORENTZ_H_

#include "core/behavior/behavior.h"
#include "core/container/math_array.h"

namespace bdm 
{
  class PenningTrapField;
  class LorentzForce : public Behavior 
  {
    BDM_BEHAVIOR_HEADER(LorentzForce, Behavior, 1);

  public:
    LorentzForce() = default;
    explicit LorentzForce(const PenningTrapField* field) : field_(field) {}
    ~LorentzForce() override = default;

    void Run(Agent* agent) override;

    const PenningTrapField* GetField() const { return field_; }
    void SetField(const PenningTrapField* field) { field_ = field; }

  private:
    const PenningTrapField* field_ = nullptr;
  };
}

#endif