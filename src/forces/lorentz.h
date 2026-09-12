#pragma once

#include "core/behavior/behavior.h"

namespace bdm 
{
    class LorentzForce : public Behavior 
    {
        BDM_BEHAVIOR_HEADER(LorentzForce, Behavior, 1);
        public:
        LorentzForce() = default;
        void Run(Agent* agent) override;
    };
} // namespace bdm