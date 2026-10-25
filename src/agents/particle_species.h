#ifndef AGENTS_PARTICLE_SPECIES_H_
#define AGENTS_PARTICLE_SPECIES_H_

#include <cstdint>
#include "utils/constants.h"

namespace bdm 
{
    enum class ParticleSpecies : uint8_t 
    {
        kAntiproton = 0,
        kPositron = 1,
        kAntihydrogen = 2,
        kElectron = 3,
        kProton = 4,
        kAntihelium = 5,
    };

    // Default physical properties (mass in kg, charge in Coulombs) associated with a given species
    struct SpeciesProperties 
    {
        real_t mass;
        real_t charge;
    };

    // Looks up the default mass/charge for a given species.
    inline SpeciesProperties GetDefaultSpeciesProperties(ParticleSpecies species) 
    {
        switch (species) 
        {
            case ParticleSpecies::kAntiproton:
            return {constants::kProtonMass, -constants::kElementaryCharge};
            case ParticleSpecies::kPositron:
            return {constants::kElectronMass, +constants::kElementaryCharge};
            case ParticleSpecies::kAntihydrogen:
            return {constants::kProtonMass + constants::kElectronMass, 0.0};
            case ParticleSpecies::kElectron:
            return {constants::kElectronMass, -constants::kElementaryCharge};
            case ParticleSpecies::kProton:
            return {constants::kProtonMass, +constants::kElementaryCharge};
            case ParticleSpecies::kAntihelium:
            return {2 * (constants::kProtonMass + constants::kNeutronMass), 0.0};
            default:
            return {0.0, 0.0};
        }
    }
}

#endif