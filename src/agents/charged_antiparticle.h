#ifndef AGENTS_CHARGED_ANTIPARTICLE_H_
#define AGENTS_CHARGED_ANTIPARTICLE_H_

#include <cstring>

#include "agents/particle_species.h"
#include "core/agent/agent.h"
#include "core/container/math_array.h"
#include "utils/constants.h"

namespace bdm 
{
  // A single charged antiparticle (or any charged particle listed in ParticleSpecies) tracked as a BioDynaMo Agent.
  class ChargedAntiparticle : public Agent 
  {
    BDM_AGENT_HEADER(ChargedAntiparticle, Agent, 4);

  public:
    ChargedAntiparticle() = default;

    // Constructs a particle of the given species at `position`, with velocity, mass, and charge initialized to the species' defaults.
    // Diameter is set to a small default for  point-charge neighbor search; can be overriden via SetDiameter if a
    // different value is needed.
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

    // Changes the agent's species and, unless `keep_overrides` is true, resets mass/charge to that species' defaults.
    void SetSpecies(ParticleSpecies species, bool keep_overrides = false);

    // --- annihilation state ---------------------------------------------------
    bool IsAnnihilated() const { return annihilated_; }
    void MarkAnnihilated() { annihilated_ = true; } // Marks the agent as annihilated

    // --- Coulomb force scratch value (see class-level comment) ---------------
    const Real3& GetCoulombForce() const { return coulomb_force_; }
    void SetCoulombForce(const Real3& force) { coulomb_force_ = force; }

    // --- Agent interface: shape, position, diameter, displacement ------------
    Shape GetShape() const override { return Shape::kSphere; }

    const Real3& GetPosition() const override { return position_; }
    void SetPosition(const Real3& pos) override { position_ = pos; }

    real_t GetDiameter() const override { return diameter_; }
    void SetDiameter(real_t diameter) override { diameter_ = diameter; }

    Real3 CalculateDisplacement(const InteractionForce* force, real_t squared_radius, real_t dt) override 
    {
      return Real3{0.0, 0.0, 0.0};
    }

    void ApplyDisplacement(const Real3& displacement) override { position_ = position_ + displacement; }

  private:
    Real3 position_ = {0.0, 0.0, 0.0};
    real_t diameter_ = 0.0;
    Real3 velocity_ = {0.0, 0.0, 0.0};
    real_t mass_ = 0.0;
    real_t charge_ = 0.0;
    ParticleSpecies species_ = ParticleSpecies::kAntiproton;
    bool annihilated_ = false;
    Real3 coulomb_force_ = {0.0, 0.0, 0.0};

  public:
    static const char* TypeName() { return "ChargedAntiparticle"; }
  };

  inline ChargedAntiparticle* FastCast(Agent* agent) 
  {
    if (agent == nullptr) { return nullptr; }
    if (std::strcmp(agent->GetTypeName(), ChargedAntiparticle::TypeName()) != 0) { return nullptr; }
    return static_cast<ChargedAntiparticle*>(agent);
  }

  inline const ChargedAntiparticle* FastCast(const Agent* agent) {
    if (agent == nullptr) { return nullptr; }
    if (std::strcmp(agent->GetTypeName(), ChargedAntiparticle::TypeName()) != 0) { return nullptr; }
    return static_cast<const ChargedAntiparticle*>(agent);
  }
}
#endif