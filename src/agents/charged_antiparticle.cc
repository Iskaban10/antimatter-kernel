#include "agents/charged_antiparticle.h"

namespace bdm 
{
    ChargedAntiparticle::ChargedAntiparticle(const Real3& position, ParticleSpecies species): species_(species) 
    {
        SetPosition(position);
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
}  // namespace bdm