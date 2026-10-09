# [Package Name]

**A Python toolkit for analysing, processing, and visualising astrophysical simulations.**

`[Package Name]` is a Python package designed to simplify the analysis of cosmological hydrodynamical simulations, with a focus on particle-based data, galaxy properties, sink particles, and black hole evolution.

The package provides a structured interface for accessing simulation outputs, constructing particle and galaxy catalogues, tracking particle evolution across snapshots, analysing black hole growth and merger histories, and producing scientific visualisations.

The goal is to make simulation analysis more modular, reusable, and reproducible by separating data access, scientific objects, catalogue construction, analysis, and visualisation.

## Features

- **Simulation management** — Centralised access to simulation metadata, paths, readers, writers, and output files.
- **Snapshot analysis** — Work with individual simulation outputs through a snapshot-based interface.
- **Particle catalogues** — Access and analyse gas cells, dark matter particles, stars, and sink particles.
- **Metadata management** — Store and retrieve particle properties, including IDs, initial masses, particle types, formation times, and stellar lifetimes.
- **Particle histories** — Track particle masses, positions, velocities, and other time-dependent properties across simulation outputs.
- **Black hole analysis** — Analyse black hole growth, accretion histories, and merger events.
- **Galaxy catalogues** — Construct and work with galaxy-level data and associated particle populations.
- **Binary data processing** — Read sink-particle information and build persistent metadata and history databases.
- **Incremental processing** — Reuse previously generated databases and process new simulation outputs without rebuilding everything from scratch.
- **Restart-aware processing** — Handle simulation restarts by identifying overlapping output times and rebuilding affected histories.
- **Scientific visualisation** — Generate projections, slices, particle renderings, and custom scientific plots using tools such as `yt` and `Matplotlib`.

## Design philosophy

The package is organised around a few core principles.

1. **Separation of concerns.** Reading raw simulation data, representing scientific objects, building catalogues, performing analysis, and creating visualisations are handled by distinct components.
2. **Lazy access where appropriate.** Expensive data structures, such as particle catalogues, can be constructed when needed rather than loaded eagerly.
3. **Reusable scientific objects.** Simulation, snapshot, particle, galaxy, and history objects provide consistent interfaces across analysis workflows.
4. **Persistent analysis products.** Expensive preprocessing steps can be saved and reused across sessions.
5. **Explicit data handling.** Missing metadata and invalid values should be handled explicitly rather than silently replaced with inferred values.

## Architecture

The package separates simulation access, scientific data structures, processing, and visualisation.

```text
[Package Name]/
├── analysis/
│   ├── core/
│   │   └── Simulation
│   │
│   ├── objects/
│   │   ├── Snapshot
│   │   ├── ParticleCatalog
│   │   ├── GalaxyCatalog
│   │   ├── Particle
│   │   ├── ParticleHistory
│   │   └── Metadata
│   │
│   ├── builders/
│   │   ├── ParticleBuilder
│   │   └── GalaxyBuilder
│   │
│   ├── inout/
│   │   ├── BinaryReader
│   │   ├── Readers
│   │   └── Writers
│   │
│   ├── visualizations/
│   │   ├── Plotter
│   │   ├── yt utilities
│   │   └── Particle renderers
│   │
│   └── ...
│
├── data/
├── examples/
├── tests/
├── README.md
└── pyproject.toml
```

*The directory tree is illustrative; update it to match the final repository layout.*

### Core components

#### `Simulation`

Represents a simulation and provides access to its configuration, directory structure, readers, writers, and available outputs.

It acts as the main entry point for working with simulation data.

#### `Snapshot`

Represents a single simulation output. It provides access to snapshot-specific information and associated scientific objects, including particle and galaxy catalogues.

Expensive objects, such as particle catalogues, can be initialised lazily when accessed.

#### `ParticleCatalog`

Provides an interface for particle populations within a snapshot. It supports working with particle properties and spatial selections without requiring each analysis script to implement its own data-access logic.

#### `GalaxyCatalog`

Provides access to galaxy-level objects and their associated properties, enabling analyses at the galaxy scale rather than directly manipulating raw particle arrays.

#### Metadata and particle histories

Metadata and histories serve different purposes:

- **Metadata** stores particle properties that are independent of a particular snapshot, such as particle ID, initial mass, formation time, and particle type.
- **Particle histories** store the evolution of individual particles, including their masses, positions, velocities, and recorded dump times.

This separation allows particle properties and time-dependent evolution to be accessed independently.

#### Builders

Builders construct persistent scientific data products from raw simulation outputs.

For example, the particle builder processes sink-particle binary dumps to create metadata databases and particle histories. It supports incremental processing and restart handling, reducing unnecessary work when additional simulation outputs become available.

#### Readers and writers

The input/output layer isolates file-format details from the rest of the package.

Readers handle simulation outputs, binary sink-particle files, and previously saved analysis products. Writers persist derived catalogues, metadata, histories, and build state for subsequent use.

#### Visualisation tools

The visualisation layer provides reusable utilities for scientific plotting and rendering, including `yt` projections and slices, Matplotlib figures, and custom particle renderers.

It is intended to support both exploratory analysis and publication-quality figures.

## Installation

Installation instructions will be added once the package's dependency configuration and distribution method are finalised.

For development, clone the repository:

```bash
git clone <repository-url>
cd <repository-directory>
```

Install the package and its dependencies using the instructions provided by the project.

**Requirements:** Python and the scientific Python ecosystem. Depending on the functionality used, dependencies may include NumPy, Matplotlib, `yt`, and the libraries required to read and process the relevant simulation output formats.

## Usage

The examples below illustrate the intended workflow. Adapt the import paths and initialisation arguments to the current public API.

### 1. Open a simulation

```python
from <package_name> import Simulation

sim = Simulation(...)

print(sim)
```

The simulation object provides the entry point for accessing the simulation's outputs and analysis tools.

### 2. Access a snapshot

```python
snapshot = sim.snapshot(...)

particles = snapshot.particles
galaxies = snapshot.galaxies
```

The snapshot provides a common interface for working with particle and galaxy data at a particular simulation time.

### 3. Analyse particle populations

```python
# Example workflow; adapt to the actual catalogue API.
stars = snapshot.particles.stars
black_holes = snapshot.particles.sinks.black_holes

print(len(black_holes))
```

Particle catalogues provide a convenient starting point for spatial selections, property analysis, and comparisons between particle populations.

### 4. Build particle metadata and histories

```python
# Example workflow; adapt to the actual builder API.
builder = ParticleBuilder(sim)

metadata, histories = builder.build()
```

The builder processes sink-particle output files and constructs reusable metadata and history databases. Existing results can be reused when no new files need processing.

### 5. Visualise simulation data

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

# Add simulation-derived data to the axes.
ax.set_xlabel("x")
ax.set_ylabel("y")

plt.show()
```

The package's visualisation utilities can be used alongside standard scientific Python tools to create projections, particle renderings, and publication figures.

## Scientific applications

The package is being developed to support analyses of cosmological simulations, particularly the formation and evolution of galaxies and black holes in the early Universe.

Potential applications include:

- Black hole mass growth and accretion histories.
- Black hole merger catalogues and merger trees.
- Stellar population and galaxy structure analysis.
- Gas density, temperature, metallicity, and velocity fields.
- Comparisons between simulations with different physical prescriptions.
- Preparation of scientific figures and derived datasets.

## Data and reproducibility

The package distinguishes between raw simulation outputs and derived analysis products.

Derived products may include particle metadata databases, particle histories, galaxy catalogues, and build-state information. Persistent outputs allow expensive preprocessing to be reused across analysis sessions.

For reliable results, users should retain the provenance of their simulation outputs and document relevant physical assumptions, units, cosmological parameters, and processing choices.

## Development status

The package is under active development. Interfaces, module organisation, and supported workflows may change as the implementation evolves.

Contributions, bug reports, and suggestions are welcome. Please include a minimal reproducible example and relevant error messages when reporting issues.

## Contributing

Contributions are welcome, particularly in the areas of:

- Testing and validation of scientific calculations.
- Support for additional simulation output formats.
- Performance and memory optimisation.
- Documentation and examples.
- Improvements to catalogue construction and visualisation.

Before submitting a contribution, please describe the intended change and ensure that existing tests pass.

## Citation

If you use this package in scientific work, please cite the repository and any associated publication once the package has an established citation record.

A `CITATION.cff` file will be added when the preferred citation information is finalised.

## License

The licensing terms will be specified in the repository's `LICENSE` file.

## Acknowledgements

This project builds on the open-source scientific Python ecosystem. Specific dependencies and relevant publications will be acknowledged in the final documentation.

