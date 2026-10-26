#include "agents/charged_antiparticle.h"

namespace bdm 
{
    ChargedAntiparticle::ChargedAntiparticle(const Real3& position, ParticleSpecies species, const Real3& velocity): species_(species)
    {
        SetPosition(position);
        SetVelocity(velocity);
        SetDiameter(constants::kDefaultParticleDiameter);

        const SpeciesProperties props = GetDefaultSpeciesProperties(species_);
        mass_ = props.mass;
        charge_ = props.charge;
    }

    void ChargedAntiparticle::SetSpecies(ParticleSpecies species, bool keep_overrides) 
    {
        species_ = species;
        if (!keep_overrides) 
        {
            const SpeciesProperties props = GetDefaultSpeciesProperties(species_);
            mass_ = props.mass;
            charge_ = props.charge;
        }
    }
}