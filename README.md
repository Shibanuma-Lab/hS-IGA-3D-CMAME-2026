# S-IGA Circular Crack in 3D Solid (Linear Analysis)

[![Version](https://img.shields.io/badge/version-1.0-blue.svg)](https://github.com/yourusername/S-IGA-circular-crack-in-3D-solid-linear)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Fortran](https://img.shields.io/badge/Fortran-90+-green.svg)](https://fortran-lang.org/)
[![Lab](https://img.shields.io/badge/Lab-Shibanuma%20Lab-red.svg)](http://www.struct.t.u-tokyo.ac.jp/shibanuma/)
[![University](https://img.shields.io/badge/University-UTokyo-orange.svg)](https://www.u-tokyo.ac.jp/)

## 📖 Overview

This repository implements a **S-version Isogeometric Analysis (S-IGA)** framework for simulating three-dimensional circular crack propagation in elastic solids. The method combines the advantages of NURBS-based Isogeometric Analysis (IGA) with local mesh refinement strategies to accurately capture crack-tip stress singularities and crack front behavior.

The implementation is a Python-based mesh generation and pre-processing framework coupled with a Fortran-based finite element solver (`sfem_linear`).

## 📁 Repository Structure

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

## 🔧 Installation

### ✅ Prerequisites

- Python 3.8+
- NumPy, SciPy
- Fortran compiler (gfortran or Intel Fortran)
- OpenMP support
- Git

### 🚀 Setup

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

## ⚡ Quick Start

### 🎯 Basic Workflow

```bash
cd circular_crack

# 1. Run simulation for steps 0-10
python3 main.py --step_start 0 --step_end 10

# 2. Calculate J-integral and stress intensity factors from results
python3 main.py --is_K --step_start 1 --step_end 10
```

### 💡 Example Use Cases

#### Example 1: 🌐 Generate Mesh Only (No Solver)
```bash
# Generate mesh and input files for step 5 without running solver
python3 main.py --step_start 5 --step_end 6 --meshonly
```

#### Example 2: ⚙️ Run Solver Only (Mesh Already Exists)
```bash
# Run solver for existing mesh files
python3 main.py --step_start 5 --step_end 6 --solveronly
```

#### Example 3: 🎲 Single Step with Fresh Start
```bash
# Run single step 20 with no restart from previous step
python3 main.py --particular --step_start 20 --no_restart
```

#### Example 4: 📊 Post-Process Results (J-Integral Calculation)
```bash
# Calculate J-integral and K_I for steps 50-100
python3 main.py --is_K --step_start 50 --step_end 100

# With custom J-integral domain parameters
python3 main.py --is_K --step_start 50 --step_end 100 \
    --Rj0 1.5 --Rj1 2.0 --Wj0 1.0 --Wj1 1.5 \
    --velocity 1200.0 --output_K results/custom_J.csv
```

#### Example 5: 🧹 Clean and Restart
```bash
# Delete all previous input files and start fresh
python3 main.py --delete --step_start 0 --step_end 10
```

## 📚 Usage

### ⚙️ Configuration

Key parameters are configured in the `circular_crack/const/` directory:

#### 🔬 Material Properties (`material_property.py`)
```python
EE = 2.06e11          # Young's modulus [Pa]
Nu = 0.3              # Poisson's ratio
Rho = 7800.0          # Density [kg/m³]
SigmaInfinity = 1.0e11  # Applied stress [Pa]
```

#### 🎮 Simulation Parameters (`simulation_params.py`)
```python
c = 4.0e-3            # Initial crack radius [m]
V = 1000.0            # Crack velocity [m/s]
step_start = 0        # Starting step
step_end = 10         # Ending step
nbcebc = 1            # BC type (0: force, 1: displacement)
```

#### 🌐 Mesh Parameters
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

- **📐 J-integral parameters** (`const_jintegral.py`):
  - Inner radius: `Rj0 = 1.5`
  - Outer radius: `Rj1 = 1.515`
  - Inner width: `Wj0 = 1.0`
  - Outer width: `Wj1 = 1.01`
  - Step range: `step_start = 1, step_end = 100`
  - Output file: `output_file = "J_integral_results.csv"`

### 💻 Command-Line Interface

The `main.py` script provides a comprehensive command-line interface for all simulation tasks:

#### 🎛️ Simulation Control
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

#### 🔬 J-Integral and Fracture Analysis
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

## 🤝 Contributing

This is a research code developed for academic purposes. Contributions should maintain:
- Code clarity and documentation
- Compatibility with existing mesh format
- Validation against analytical solutions

## 📚 References

1. Hughes, T.J.R. et al. "Isogeometric Analysis: CAD, Finite Elements, NURBS, Exact Geometry and Mesh Refinement." *Computer Methods in Applied Mechanics and Engineering*, 2005.

## 📄 License

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

## 👨‍💻 Authors

**Tianyu He**
- Primary Developer and Maintainer

## 🙏 Acknowledgments

This work uses the `sfem_linear` solver framework and builds upon established methods in computational fracture mechanics and isogeometric analysis.
