# Antimatter-Kernel

An agent-based simulation of charged antiparticles confined in a Penning trap, built on BioDynaMo.

The project follows the [Agent-based Large-Scale Antimatter Simulation](https://compiler-research.org/open_projects?tag=all#agent-based-large-scale-antimatter-simulation) project from compiler-research.org, which aims to extend BioDynaMo to charged-particle simulation of the kind used in antimatter experiments such as CERN's AEgIS.

>Status: Early Development


## Table of Contents
1. [Building the simulation](#building-the-simulation)
2. [Running the simulation](#running-the-simulation)
3. [Project Overview](#project-overview)
4. [Project Structure](#project-structure)
5. [Input Parameters](#input-parameters)

---
## Building the simulation

Option 1:
```bash
biodynamo build
```

Option 2:
```bash
mkdir build && cd build
cmake ..
make -j <number_of_processes_for_build>
```

## Running the simulation

Before running the simulation, you need to build it. If you haven't done so, 
please go back to step 1.

Option 1:
```bash
biodynamo run
```

Option 2:
```
./build/antimatter-kernel
```
---
--- 

## Project Overview

AntimatterKernel is an agent-based model, written in C++ on top of BioDynaMo, that simulates such ensembles. Each antiparticle is an agent carrying its own position, velocity, mass, charge and species, and the module provides the physics as composable BioDynaMo plugins:

- **Penning-trap external fields:** a uniform axial magnetic field plus a hyperbolic-electrode quadrupole electrostatic potential.
- **Lorentz-force dynamics:** particle motion is advanced with the Boris algorithm, the standard integrator for charged-particle motion in magnetic fields. Unlike explicit Euler integration, it does not introduce spurious numerical heating of the cyclotron motion over many orbits.
- **Pairwise Coulomb interaction:** computed through BioDynaMo's neighbor-search infrastructure with a configurable cutoff radius. The Coulomb force is folded into the same symmetric half-kick as the trap's electric field inside the Boris step, so the combined system keeps second-order accuracy instead of being applied as a separate, less accurate velocity kick.
- **Stochastic annihilation:** particle-antiparticle pairs (positron-electron, antiproton-proton) that come into contact can annihilate with a configurable probability, using a two-step mark-then-remove scheme so annihilation events remain observable in the output for one step.
- **Multi-species support:** antiproton, positron, antihydrogen, electron and proton, each with SI-unit default mass and charge.

The project aims to turn Penning-trap prototypes into a reusable, documented and scalable module for antimatter-motivated studies and other charged-particle systems.

### Validation: Energy Analysis

The magnetic field does no work on a charged particle, so for a correct and adequately time-resolved integrator the total energy (kinetic + trap potential + Coulomb potential) should stay approximately constant during a run. The module logs this at regular intervals and writes it to `energy_analysis.csv` inside the run's output directory. Steady, monotonic growth of the fractional drift is the signature of an integration problem (for example `dt` too large relative to the cyclotron period). The Coulomb potential energy is computed with the *same* cutoff radius as the force, so any drift observed reflects the integrator and not the cutoff approximation.

---

## Project Structure

The project is organized into the following components:

### Core Simulation Files (`src/`)

- **[`antimatter-kernel.h`](src/antimatter-kernel.h) / [`antimatter-kernel.cc`](src/antimatter-kernel.cc)**: The entry point. The header contains the inline `Simulate` function, which registers the parameter groups, loads `params.json`, builds the Penning-trap field, seeds the particle cloud, attaches the behaviors to each agent, runs the simulation step by step, and logs the energy analysis. The `.cc` file only contains `main`.

- **[`charged_antiparticle.h`](src/agents/charged_antiparticle.h) / [`charged_antiparticle.cc`](src/agents/charged_antiparticle.cc)**: The `ChargedAntiparticle` agent. It stores position, diameter, velocity, mass, charge, species and annihilation state, implements the virtual position, diameter, displacement interface required by BioDynaMo's `Agent`, and holds a per-step Coulomb-force value consumed by the Lorentz force.

- **[`particle_species.h`](src/agents/particle_species.h)**: The `ParticleSpecies` enum and a table of default SI mass and charge per species.

- **[`penning_trap.h`](src/fields/penning_trap.h) / [`penning_trap.cc`](src/fields/penning_trap.cc)**: `PenningTrapField`, which returns the electric field, magnetic field and electric potential at a given position. It consists of pure field geometry i.e. it knows nothing about any particle's mass, charge or velocity.

- **[`lorentz.h`](src/forces/lorentz.h) / [`lorentz.cc`](src/forces/lorentz.cc)**: The `LorentzForce` behavior. Advances velocity and position with the Boris algorithm, combining the trap electric field and the stored Coulomb force in one symmetric half-kick around the magnetic rotation. It reads `dt` from the registered `PhysicsParams`.

- **[`coulomb_interaction.h`](src/interaction/coulomb_interaction.h) / [`coulomb_interaction.cc`](src/interaction/coulomb_interaction.cc)**: The `CoulombInteraction` behavior. Sums the pairwise Coulomb force on each agent from neighbors within the cutoff radius (via BioDynaMo's `Environment::ForEachNeighbor`) and stores it on the agent. It does not integrate velocity or position itself. **It must run before `LorentzForce` on each agent.**

- **[`annihilation.h`](src/interaction/annihilation.h) / [`annihilation.cc`](src/interaction/annihilation.cc)**: The `AnnihilationInteraction` behavior and the `IsAntiparticlePair` helper. Marks matching particle-antiparticle pairs as annihilated with the configured probability, then removes them from the simulation on the following step.

- **[`physics_params.h`](src/params/physics_params.h) / [`physics_params.cc`](src/params/physics_params.cc)**, **[`interaction_params.h`](src/params/interaction_params.h) / [`interaction_params.cc`](src/params/interaction_params.cc)**, **[`stochastic_params.h`](src/params/stochastic_params.h) / [`stochastic_params.cc`](src/params/stochastic_params.cc)**: These three BioDynaMo `ParamGroup`s hold the trap/time-integration/initial-cloud settings, the Coulomb settings, and the annihilation settings respectively. Each provides `LoadParams` (read from `params.json`, keeping defaults for anything not specified), `PrintParams`, and `Validate`.

- **[`energy_analysis.h`](src/analysis/energy_analysis.h) / [`energy_analysis.cc`](src/analysis/energy_analysis.cc)**: It contains `ComputeTotalEnergy`, which returns the kinetic, trap-potential and Coulomb-potential energy and the number of active (non-annihilated) particles.

- **[`constants.h`](src/utils/constants.h)**: Physical constants in SI units and shared simulation defaults, typed as BioDynaMo's `real_t`.

### Configuration Files

- **[`params.json`](params.json)**: JSON file for configuring physics, interaction and annihilation settings without recompiling.

- **[`bdm.toml`](bdm.toml)**: BioDynaMo-specific configuration for the neighbor-search backend, output directory and visualization export.

### Output

Each run writes into BioDynaMo's per-run output directory (`output/antimatter-kernel/`):

- **`energy_analysis.csv`**: one row per logged step with the columns `step`, `kinetic_energy_J`, `trap_potential_energy_J`, `coulomb_potential_energy_J`, `total_energy_J`, `num_active_particles` and `fractional_drift` (the relative change of total energy from step 0).
---

## Input Parameters

All simulation parameters are set in `params.json` without the need to recompile. Any parameter that is not specified keeps its default value. Parameters are grouped by the `ParamGroup` that reads them.

### General Simulation Parameters (`PhysicsParams`)

These are the parameters most commonly changed when designing a run.

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `seed` | 42 | - | Seed for the random placement of the initial particle cloud, for reproducibility |
| `species` | `"positron"` | - | Species of every seeded particle: `antiproton`, `positron`, `antihydrogen`, `electron` or `proton` |
| `num_particles` | 10 | - | Number of particles seeded at the start |
| `initial_cloud_radius` | 1e-4 | m | Particles are placed at random positions within ± this distance of the trap center. Must be well inside `z0` and `rho0` |
| `dt` | 1e-12 | s | Integration timestep. Must be much smaller than the cyclotron period 2πm / (\|q\|B₀) (about 36 steps per period by default) |
| `num_steps` | 50000 | - | Number of simulation steps. Total simulated time is `dt × num_steps` |

### Trap Parameters (`PhysicsParams`)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `b0` | 1.0 | T | Strength of the uniform axial magnetic field |
| `trap_voltage` | 10.0 | V | Trap (ring-to-endcap) voltage U₀ of the quadrupole potential. A positive value axially confines a positively charged species; use a negative value to trap a negative species |
| `z0` | 0.01 | m | Endcap half-distance |
| `rho0` | 0.01 | m | Ring-electrode radius |

The characteristic trap dimension is derived internally as `d² = z0²/2 + rho0²/4`.

### Coulomb Interaction Parameters (`InteractionParams`)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `coulomb_enabled` | true | - | Enables or disables the pairwise Coulomb interaction. Disable it to isolate single-particle (Lorentz-only) validation runs |
| `cutoff_radius` | 1e-3 | m | Agents farther apart than this do not exert Coulomb force on each other. **There is no built-in default in the code (0 is rejected as invalid when Coulomb is enabled), so a value must be provided.** The shipped value is chosen to be much larger than the expected inter-particle spacing and much smaller than the trap dimensions; it should be revisited for different densities |
| `barnes_hut_theta` | 0.5 | - | Sets the theta for the Barnes-Hut algorithm implementation, lies between 0 and 1 |
| `fast_mulitpole_order` | 0.5 | - | Sets the order for the Fast-Multipole algorithm implementation, higher the value more accurate the results |
| `fast_multipole_theta` | 0.5 | - | Sets the theta for the Fast-Multipole algorithm implementation, lies between 0 and 1 |

### Annihilation Parameters (`StochasticParams`)

| Parameter | Default | Units | Description |
|-----------|---------|-------|-------------|
| `annihilation_enabled` | false | - | Enables or disables stochastic annihilation. Keep disabled unless specifically studying the mechanism |
| `contact_radius` | 1e-9 | m | A particle-antiparticle pair closer than this is considered in contact and an annihilation attempt is evaluated |
| `annihilation_probability` | 0.1 | - | Probability in [0, 1] of annihilation per step, per pair found within `contact_radius` |

### Parameter Validation

At startup every parameter group is validated and the program stops with an error message on an invalid configuration, for example `b0 <= 0`, `dt <= 0`, `initial_cloud_radius` not smaller than `z0` and `rho0`, `cutoff_radius <= 0` while Coulomb is enabled, or `annihilation_probability` outside [0, 1]. An unrecognized `species` name also stops the run rather than silently falling back to a default. The effective values used are printed at the start of every run.

### Parameter Configuration

To modify parameters, edit [`params.json`](params.json). Only include the parameters you want to change from their defaults. Example:

```json
{
  "seed": 42,
  "species": "positron",
  "num_particles": 10,
  "dt": 1e-12,
  "num_steps": 50000,
  "b0": 1.0,
  "trap_voltage": 10.0,

  "coulomb_enabled": true,
  "coulomb_method": "standard",
  "cutoff_radius": 1e-3,
  "barnes_hut_theta": 0.5,
  "fast_multipole_order": 4,
  "fast_multipole_theta": 0.5,

  "annihilation_enabled": false
}
```

### BioDynaMo Engine Configuration (`bdm.toml`)

| Key | Example | Description |
|-----|---------|-------------|
| `[simulation] output_dir` | `"output"` | Base directory for simulation output |
| `[simulation] environment` | `"kd_tree"` | Neighbor-search backend: `uniform_grid`, `kd_tree` or `octree` |
| `[visualization] export` | `true` | Export particle data for ParaView |
| `[visualization] export_generate_pvsm` | `true` | Generate a ParaView state file (`.pvsm`) automatically |
| `[visualization] interval` | `10` | Export every N-th step. With `num_steps = 50000` and `interval = 1` this writes 50,000 file sets, so raise it for long runs |
| `[[visualize_agent]] name` | `"ChargedAntiparticle"` | Agent type to export (one block per agent type) |