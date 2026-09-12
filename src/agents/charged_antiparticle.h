#ifndef AGENTS_CHARGED_ANTIPARTICLE_H_
#define AGENTS_CHARGED_ANTIPARTICLE_H_

#include "agents/particle_species.h"
#include "core/agent/agent.h"
#include "core/container/math_array.h"
#include "utils/constants.h"

namespace bdm 
{
  class ChargedAntiparticle : public Agent 
  {
      BDM_AGENT_HEADER(ChargedAntiparticle, Agent, 3);

    public:
      ChargedAntiparticle() = default;

      /// Constructs a particle at `position`, with mass and charge initialized to the species' defaults with 0 velocity
      explicit ChargedAntiparticle(const Real3& position, ParticleSpecies species = ParticleSpecies::kAntiproton);

      ~ChargedAntiparticle() override = default;

      // --- kinematics ---------------------------------------------------------
      const Real3& GetVelocity() const { return velocity_; }
      void SetVelocity(const Real3& velocity) { velocity_ = velocity; }

      // --- physical properties -------------------------------------------------
      real_t GetMass() const { return mass_; }
      void SetMass(real_t mass) { mass_ = mass; }

      real_t GetCharge() const { return charge_; }
      void SetCharge(real_t charge) { charge_ = charge; }

      ParticleSpecies GetSpecies() const { return species_; }
      void SetSpecies(ParticleSpecies species, bool keep_overrides = false);

      // --- annihilation state ---------------------------------------------------
      bool IsAnnihilated() const { return annihilated_; }
      void MarkAnnihilated() { annihilated_ = true; } // marks the agent as annihilated.

      // --- Agent interface: shape, position, diameter, displacement ------------
      Shape GetShape() const override { return Shape::kSphere; }

      const Real3& GetPosition() const override { return position_; }
      void SetPosition(const Real3& pos) override { position_ = pos; }

      real_t GetDiameter() const override { return diameter_; }
      void SetDiameter(real_t diameter) override { diameter_ = diameter; }

      // override bdm's calculate displacement with default return values
      Real3 CalculateDisplacement(const InteractionForce* force, real_t squared_radius, real_t dt) override { return Real3{0.0, 0.0, 0.0}; }
      void ApplyDisplacement(const Real3& displacement) override { position_ = position_ + displacement; }

    private:
      Real3 position_ = {0.0, 0.0, 0.0};
      real_t diameter_ = 0.0;
      Real3 velocity_ = {0.0, 0.0, 0.0};
      real_t mass_ = 0.0;
      real_t charge_ = 0.0;
      ParticleSpecies species_ = ParticleSpecies::kAntiproton;
      bool annihilated_ = false;
  };
}  // namespace bdm
#endif  // AGENTS_CHARGED_ANTIPARTICLE_H_