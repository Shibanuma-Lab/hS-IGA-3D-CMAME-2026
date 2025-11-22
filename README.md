# S-IGA Circular Crack in 3D Solid (Linear Analysis)

## Overview

This repository implements a **S-version Isogeometric Analysis (S-IGA)** framework for simulating three-dimensional circular crack propagation in elastic solids. The method combines the advantages of NURBS-based Isogeometric Analysis (IGA) with local mesh refinement strategies to accurately capture crack-tip stress singularities and crack front behavior.

The implementation is a Python-based mesh generation and pre-processing framework coupled with a Fortran-based finite element solver (`sfem_linear`).

## Repository Structure

```
S-IGA-circular-crack-in-3D-solid-linear/
├── circular_crack/              # Main Python package
│   ├── main.py                  # Entry point with command-line interface
│   ├── global_mesh.py           # NURBS-based global mesh generation
│   ├── local_mesh.py            # Local mesh refinement around crack
│   ├── boundary.py              # Boundary condition generation (EBC/NBC)
│   ├── initial.py               # Initial condition setup for restart
│   ├── input_generator.py       # Input file generation for solver
│   ├── jintegral.py             # J-integral and stress intensity factor calculation
│   ├── const/                   # Configuration parameters
│   │   ├── simulation_params.py # Simulation control parameters
│   │   ├── material_property.py # Material properties (E, ν, ρ, etc.)
│   │   ├── const_global_mesh.py # Global mesh parameters
│   │   ├── const_local_mesh.py  # Local mesh parameters
│   │   └── const_jintegral.py   # J-integral calculation parameters
│   ├── utils/                   # Utility functions
│   │   ├── logger.py            # Logging utility
│   │   └── step2str.py          # Step number formatting
│   └── scripts/                 # Execution scripts
│       └── linux_command.py     # Solver execution interface
├── sfem_linear/                 # Fortran solver
│   ├── bin/sfem_linear          # Compiled solver executable
│   ├── src/                     # Fortran source code
│   ├── example/                 # Example cases
│   └── manual/                  # Solver documentation
└── README.md                    # This file
```

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

## Quick Start

### Basic Workflow

```bash
cd circular_crack

# 1. Run simulation for steps 0-10
python3 main.py --step_start 0 --step_end 10

# 2. Calculate J-integral and stress intensity factors from results
python3 main.py --is_K --step_start 1 --step_end 10
```

### Example Use Cases

#### Example 1: Generate Mesh Only (No Solver)
```bash
# Generate mesh and input files for step 5 without running solver
python3 main.py --step_start 5 --step_end 6 --meshonly
```

#### Example 2: Run Solver Only (Mesh Already Exists)
```bash
# Run solver for existing mesh files
python3 main.py --step_start 5 --step_end 6 --solveronly
```

#### Example 3: Single Step with Fresh Start
```bash
# Run single step 20 with no restart from previous step
python3 main.py --particular --step_start 20 --no_restart
```

#### Example 4: Post-Process Results (J-Integral Calculation)
```bash
# Calculate J-integral and K_I for steps 50-100
python3 main.py --is_K --step_start 50 --step_end 100

# With custom J-integral domain parameters
python3 main.py --is_K --step_start 50 --step_end 100 \
    --Rj0 1.5 --Rj1 2.0 --Wj0 1.0 --Wj1 1.5 \
    --velocity 1200.0 --output_K results/custom_J.csv
```

#### Example 5: Clean and Restart
```bash
# Delete all previous input files and start fresh
python3 main.py --delete --step_start 0 --step_end 10
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
  - Crack surface elements: `aL = 51`
  - Ligament elements: `lL = 15`
  - Thickness elements: `HL = 11`
  - Angular resolution: `d_theta = 3°`

- **J-integral parameters** (`const_jintegral.py`):
  - Inner radius: `Rj0 = 1.5`
  - Outer radius: `Rj1 = 1.515`
  - Inner width: `Wj0 = 1.0`
  - Outer width: `Wj1 = 1.01`
  - Step range: `step_start = 1, step_end = 100`
  - Output file: `output_file = "J_integral_results.csv"`

### Command-Line Interface

The `main.py` script provides a comprehensive command-line interface for all simulation tasks:

#### Simulation Control
```bash
python3 main.py [OPTIONS]

Options:
  --step_start N        Starting step number (default: 0)
  --step_end N          Ending step number (default: 101)
  --meshonly            Generate mesh only, skip solver
  --solveronly          Run solver only (mesh must exist)
  --particular          Run single step (step_start only)
  --restart 0/1         Force restart value (0: fresh, 1: restart)
  --no_restart          Force fresh start for all steps
  --delete              Delete inputfiles/ before running
  --debugmode           Enable verbose debug logging
```

#### J-Integral and Fracture Analysis
```bash
python3 main.py --is_K [OPTIONS]

J-Integral Options:
  --is_K                Enable J-integral calculation mode
  --Rj0 FLOAT           Inner radius for J-domain (default: 1.5)
  --Rj1 FLOAT           Outer radius for J-domain (default: 1.515)
  --Wj0 FLOAT           Inner width parameter (default: 1.0)
  --Wj1 FLOAT           Outer width parameter (default: 1.01)
  --velocity FLOAT      Crack velocity in m/s (default: 1000.0)
  --output_K FILE       Output CSV file (default: J_integral_results.csv)
```

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

#### Multi-Step Simulation via Command Line
```bash
# Run steps 0-10 using command line
python3 main.py --step_start 0 --step_end 10
```

This will:
1. Generate global and local meshes for each step
2. Apply boundary conditions
3. Create input files for the solver
4. Execute the solver
5. Save results to `circular_crack/inputfiles/step#####/`

### J-Integral and Stress Intensity Factor Calculation

After running the simulation, calculate fracture mechanics parameters:

```bash
# Basic J-integral calculation
python3 main.py --is_K --step_start 1 --step_end 100
```

This will:
1. Load displacement, velocity, and acceleration data from `results/step#####/log/`
2. Generate full-circle mesh (0-180°) from half-circle data (0-90°)
3. Calculate J-integral for each crack front angle (0° to 90° in 3° increments)
4. Compute stress intensity factors K_I using dynamic crack theory
5. Output results to CSV files:
   - `J_integral_results.csv`: J-integral values for all steps and angles
   - `J_integral_results_KId.csv`: Stress intensity factors K_I

#### J-Integral Theory

The implementation follows the domain integral method for dynamic crack propagation:

**J-Integral (Dynamic)**:
$$J = J_s + J_d = \int_{\Gamma} \left[ W \delta_{1j} - \sigma_{ij} \frac{\partial u_i}{\partial x_1} \right] q_{,j} \, dV + \int_{\Gamma} \rho \ddot{u}_i \frac{\partial u_i}{\partial x_1} q \, dV$$

where:
- $J_s$: Static component (strain energy release rate)
- $J_d$: Dynamic component (kinetic energy contribution)
- $W = \frac{1}{2} \sigma_{ij} \varepsilon_{ij}$: Strain energy density
- $q$: Weight function (1 at crack tip, 0 at domain boundary)
- $\rho$: Material density
- $\ddot{u}_i$: Acceleration field

**Stress Intensity Factor**:
$$K_I = \sqrt{\frac{E J}{(1+\nu) A_I}}$$

where $A_I$ is a dynamic correction factor depending on crack velocity:
$$A_I = \frac{\beta_1 (1 - \beta_2^2)}{4\beta_1\beta_2 - (1 + \beta_2^2)^2}$$

with $\beta_1 = \sqrt{1 - (v/c_1)^2}$, $\beta_2 = \sqrt{1 - (v/c_2)^2}$, and $c_1$, $c_2$ are wave speeds.

#### Output Format

The CSV files contain:
- **Column 1**: Step number
- **Columns 2-32**: J or K_I values at angles 0°, 3°, 6°, ..., 87°, 90° (31 angles)

Example output:
```csv
step,0.0,3.0,6.0,9.0,...,87.0,90.0
1,1.556e7,1.554e7,1.548e7,...,1.540e7,1.556e7
2,2.294e6,2.294e6,2.296e6,...,2.294e6,2.294e6
```

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

### Simulation Results

Solver output is stored in `results/step#####/log/`:

| File | Description |
|------|-------------|
| `u_gl.l.dat` | Global-local displacement field |
| `v_gl.l.dat` | Velocity field |
| `a_gl.l.dat` | Acceleration field |
| `log.txt` | Solver convergence information |

### J-Integral Results

Post-processing output:

| File | Description |
|------|-------------|
| `J_integral_results.csv` | J-integral values (all steps × 31 angles) |
| `J_integral_results_KId.csv` | Stress intensity factors K_I |

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
