#pragma once
// -----------------------------------------------------------------------
// Fundamental physical constants (SI units: kg, C, m, s).
// -----------------------------------------------------------------------
namespace bdm::constants
{
    // Elementary charge, in Coulombs (C).
    inline constexpr double kElementaryCharge = 1.602176634e-19;

    // Electron (and positron) rest mass, in kilograms (kg).
    inline constexpr double kElectronMass = 9.1093837015e-31;

    // Proton (and antiproton) rest mass, in kilograms (kg).
    inline constexpr double kProtonMass = 1.67262192369e-27;

    // Vacuum permittivity, epsilon_0, in Farads per meter (F/m).
    inline constexpr double kVacuumPermittivity = 8.8541878128e-12;

    // Coulomb's constant, k_e = 1 / (4 * pi * epsilon_0), in N*m^2/C^2.
    // Precomputed to avoid repeating the division in hot interaction loops.
    inline constexpr double kCoulombConstant = 8.9875517923e9;

    // Boltzmann constant, in Joules per Kelvin (J/K).
    inline constexpr double kBoltzmannConstant = 1.380649e-23;

    // Speed of light in vacuum, in meters per second (m/s).
    inline constexpr double kSpeedOfLight = 2.99792458e8;

    // Default agent diameter, in meters (m), used to size BioDynaMo's spatial decomposition grid for neighbor search.
    inline constexpr double kDefaultParticleDiameter = 1e-6;

}  // namespace bdm::constants