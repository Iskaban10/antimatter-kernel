#ifndef UTILS_CONSTANTS_H_
#define UTILS_CONSTANTS_H_

#include "core/util/type.h"

namespace bdm 
{
    namespace constants 
    {
        // Elementary charge, in Coulombs (C).
        inline constexpr real_t kElementaryCharge = 1.602176634e-19;

        // Electron (and positron) rest mass, in kilograms (kg).
        inline constexpr real_t kElectronMass = 9.1093837015e-31;

        // Proton (and antiproton) rest mass, in kilograms (kg).
        inline constexpr real_t kProtonMass = 1.67262192369e-27;

        // Vacuum permittivity, epsilon_0, in Farads per meter (F/m).
        inline constexpr real_t kVacuumPermittivity = 8.8541878128e-12;

        // Coulomb's constant, k_e = 1 / (4 * pi * epsilon_0), in N*m^2/C^2.
        inline constexpr real_t kCoulombConstant = 8.9875517923e9;

        // Boltzmann constant, in Joules per Kelvin (J/K).
        inline constexpr real_t kBoltzmannConstant = 1.380649e-23;

        // Speed of light in vacuum, in meters per second (m/s).
        inline constexpr real_t kSpeedOfLight = 2.99792458e8;

        // Default agent diameter, in meters (m)
        inline constexpr real_t kDefaultParticleDiameter = 1e-6;
        
    }  // namespace constants
}  // namespace bdm

#endif  // UTILS_CONSTANTS_H_