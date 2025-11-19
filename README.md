# S-IGA Circular Crack in 3D Solid (Linear Analysis)

## Overview

This repository implements a **Scaled Isogeometric Analysis (S-IGA)** framework for simulating three-dimensional circular crack propagation in elastic solids under linear mechanical loading. The method combines the advantages of NURBS-based Isogeometric Analysis (IGA) with local mesh refinement strategies to accurately capture crack-tip stress singularities and crack front behavior.

The implementation is a Python-based mesh generation and pre-processing framework coupled with a Fortran-based finite element solver (`sfem_linear`).

## Features

- **NURBS-Based Global Mesh**: Leverages B-spline basis functions for accurate geometric representation and higher-order continuity
- **Local Mesh Refinement**: Adaptive mesh refinement around crack tips using cylindrical coordinate systems
- **Crack Propagation Simulation**: Step-by-step simulation of circular crack growth with configurable crack velocity
- **Essential and Natural Boundary Conditions**: Support for both displacement-controlled (essential) and force-controlled (natural) boundary conditions
- **FEM Data Integration**: Interpolation of boundary conditions from pre-computed FEM solutions
- **Parallel Computing Support**: OpenMP-enabled solver for efficient large-scale computations

## Repository Structure

```
S-IGA-circular-crack-in-3D-solid-linear/
├── circular_crack/              # Main Python package
│   ├── main.py                  # Entry point for mesh generation and simulation
│   ├── global_mesh.py           # NURBS-based global mesh generation
│   ├── local_mesh.py            # Local mesh refinement around crack
│   ├── boundary.py              # Boundary condition generation (EBC/NBC)
│   ├── initial.py               # Initial condition setup for restart
│   ├── input_generator.py       # Input file generation for solver
│   ├── fem_data_loader.py       # FEM data loading and interpolation
│   ├── const/                   # Configuration parameters
│   │   ├── simulation_params.py # Simulation control parameters
│   │   ├── material_property.py # Material properties
│   │   ├── const_global_mesh.py # Global mesh parameters
│   │   └── const_local_mesh.py  # Local mesh parameters
│   ├── utils/                   # Utility functions
│   │   ├── logger.py            # Logging utility
│   │   └── step2str.py          # Step number formatting
│   ├── scripts/                 # Execution scripts
│   │   └── linux_command.py     # Solver execution interface
│   └── data/                    # FEM reference data
│       └── FEMdata/             # Pre-computed FEM solutions
├── sfem_linear/                 # Fortran solver (submodule)
│   ├── bin/sfem_linear          # Compiled solver executable
│   ├── src/                     # Fortran source code
│   ├── example/                 # Example cases
│   └── manual/                  # Solver documentation
└── README.md                    # This file
```

## Theoretical Background

### Scaled Isogeometric Analysis (S-IGA)

S-IGA extends traditional IGA by introducing a scaling strategy that combines:

1. **Global NURBS Mesh**: Coarse mesh using Non-Uniform Rational B-Splines (NURBS) for the far-field region
2. **Local Refined Mesh**: Fine mesh in cylindrical coordinates around the crack tip to capture singularities
3. **Visualization Mesh**: Linear hexahedral elements for post-processing and boundary condition mapping

The method achieves high accuracy in fracture mechanics simulations while maintaining computational efficiency through adaptive refinement.

### Crack Geometry

The implementation focuses on **circular penny-shaped cracks** in 3D elastic solids:
- Crack radius: `c` (default: 4 mm)
- Crack plane: z = 0
- Crack propagation: Radially outward with velocity `V`
- Crack front discretization: Polar coordinates with angular resolution `d_theta`

## Installation

### Prerequisites

- Python 3.8+
- NumPy, SciPy
- Fortran compiler (gfortran or Intel Fortran)
- OpenMP support
- Git

### Setup

1. Clone the repository:
```bash
git clone <repository-url>
cd S-IGA-circular-crack-in-3D-solid-linear
```

2. Build the Fortran solver:
```bash
cd sfem_linear
make
cd ..
```

3. Verify the solver is compiled:
```bash
ls -lh sfem_linear/bin/sfem_linear
```

## Usage

### Configuration

Key parameters are configured in the `circular_crack/const/` directory:

#### Material Properties (`material_property.py`)
```python
EE = 2.06e11          # Young's modulus [Pa]
Nu = 0.3              # Poisson's ratio
Rho = 7800.0          # Density [kg/m³]
SigmaInfinity = 1.0e11  # Applied stress [Pa]
```

#### Simulation Parameters (`simulation_params.py`)
```python
c = 4.0e-3            # Initial crack radius [m]
V = 1000.0            # Crack velocity [m/s]
step_start = 0        # Starting step
step_end = 10         # Ending step
nbcebc = 1            # BC type (0: force, 1: displacement)
```

#### Mesh Parameters
- **Global mesh** (`const_global_mesh.py`):
  - NURBS degrees: `p=2, q=2, r=2`
  - Element size ratio: `rGL = 6`
  - Mesh scaling: `mu_G = 0.99^0.5`

- **Local mesh** (`const_local_mesh.py`):
  - Element size: `hL = 0.04 mm`
  - Crack surface elements: `aL = 9`
  - Ligament elements: `lL = 15`
  - Thickness elements: `HL = 8`
  - Angular resolution: `d_theta = 3°`

### Running Simulations

#### Single Step Execution
```python
from circular_crack.main import makemeshs
from circular_crack.scripts import linux_command

# Generate mesh for step 0
step = 0
REstart = 0  # 0: new simulation, 1: restart from previous
l, g = makemeshs(step, REstart)

# Run solver
linux_command.run(step)
```

#### Multi-Step Simulation
```python
from circular_crack.main import main

# Run full simulation (step_start to step_end)
main()
```

This will:
1. Generate global and local meshes for each step
2. Apply boundary conditions
3. Create input files for the solver
4. Execute the solver
5. Save results to `circular_crack/inputfiles/step#####/`

### Output Files

For each simulation step, the following files are generated in `inputfiles/step#####/`:

| File | Description |
|------|-------------|
| `node.g.dat` | Global mesh control point coordinates |
| `elem.g.dat` | Global mesh element connectivity |
| `node.l.dat` | Local mesh node coordinates |
| `elem.l.dat` | Local mesh element connectivity |
| `weights.g.dat` | NURBS weights for global mesh |
| `index.g.dat` | Global-local mesh mapping indices |
| `bc.g.dat` | Essential boundary conditions (global) |
| `bc.l.dat` | Natural boundary conditions (local) |
| `load.dat` | Applied load configuration |
| `input.dat` | Solver control parameters |
| `node.v.dat` | Visualization mesh nodes |
| `elem.v.dat` | Visualization mesh elements |

### Results

Solver output includes:
- Displacement fields (`.dat` files)
- Stress and strain fields
- Convergence logs in `logs/`
- VTU files for visualization (if enabled)

## Technical Details

### Mesh Generation Pipeline

1. **Global Mesh** (`global_mesh.py`):
   - Generates NURBS control points using tensor product B-splines
   - Applies graded refinement near boundaries
   - Creates visualization mesh for boundary mapping
   - Outputs: `node.g.dat`, `elem.g.dat`, `weights.g.dat`, `node.v.dat`, `elem.v.dat`

2. **Local Mesh** (`local_mesh.py`):
   - Constructs cylindrical mesh around crack front
   - Applies radial refinement towards crack tip
   - Handles crack surface and ligament discretization
   - Outputs: `node.l.dat`, `elem.l.dat`

3. **Boundary Conditions** (`boundary.py`):
   - **Essential BC (EBC)**: Displacement constraints on outer boundaries
   - **Natural BC (NBC)**: Traction/force conditions on crack surfaces
   - Uses FEM data interpolation for realistic loading
   - Outputs: `bc.g.dat`, `bc.l.dat`

4. **Global-Local Coupling** (`input_generator.py`):
   - Generates mapping between global and local meshes
   - Creates index arrays for domain decomposition
   - Outputs: `index.g.dat`

### Coordinate Systems

- **Global**: Cartesian coordinates (x, y, z) for far-field domain
- **Local**: Cylindrical coordinates (r, θ, z) centered at crack center
- **Transformation**: FEM → IGA coordinates with rotation and translation

### Boundary Condition Types

#### Essential Boundary Conditions (nbcebc = 1)
Applied on outer boundaries (x_max, y_max, z_max) using interpolated displacements from FEM solutions.

#### Natural Boundary Conditions (nbcebc = 0)
Applied as nodal forces on crack surfaces and interior boundaries.

### Crack Propagation Algorithm

For each step `n`:
1. Update crack radius: `c_n = c_0 + n * V * dt`
2. Regenerate local mesh around expanded crack
3. Update boundary conditions for new crack configuration
4. Solve equilibrium equations
5. Extract crack-tip fields (stress intensity factors, J-integral)

## Solver Integration

The Python framework generates input files for the `sfem_linear` Fortran solver:

- **Solver type**: Implicit static or dynamic analysis
- **Matrix solver**: Direct (monolis) or iterative
- **Nonlinearity**: Linear elastic (current), plasticity (future)
- **Parallelization**: OpenMP threading (24 threads default)

### Input File Format (`input.dat`)

```
solution_type     # 0: static, 1: dynamic
penalty_parameter # Large value for constraint enforcement
time_step
total_time
alpha, beta, gamma  # HHT time integration parameters
n_global_nodes
n_local_nodes
n_global_elements
n_local_elements
p, q, r            # NURBS degrees
ngp                # Gauss integration points
material_params    # E, nu, rho, sigma_y
output_options
```

## Contributing

This is a research code developed for academic purposes. Contributions should maintain:
- Code clarity and documentation
- Compatibility with existing mesh format
- Validation against analytical solutions

## References

1. Hughes, T.J.R. et al. "Isogeometric Analysis: CAD, Finite Elements, NURBS, Exact Geometry and Mesh Refinement." *Computer Methods in Applied Mechanics and Engineering*, 2005.

## License

MIT License

Copyright (c) 2025 Tianyu He

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Authors

**Tianyu He**
- Primary Developer and Maintainer

## Acknowledgments

This work uses the `sfem_linear` solver framework and builds upon established methods in computational fracture mechanics and isogeometric analysis.
