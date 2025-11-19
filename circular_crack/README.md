# Circular Crack in 3D Solid - Python Implementation

This Python framework generates input files for S-IGA (Smoothed Isogeometric Analysis) simulation of circular crack propagation in a 3D solid under dynamic loading.

## Overview

This implementation is converted from Mathematica notebooks and creates:
- Global mesh with graded refinement
- Local mesh around crack tip
- Boundary conditions and loads
- Input files for SFEM linear solver
- Restart capabilities for dynamic crack propagation

## Directory Structure

```
circular_crack/
├── main.py                 # Main execution script
├── global_mesh.py          # Global mesh generation
├── local_mesh.py           # Local mesh (crack region)
├── boundary.py             # Boundary conditions and loads
├── input_generator.py      # SFEM input file generation
├── initial.py              # Initial conditions for restart
├── linux_command.py        # Solver execution
├── const/
│   ├── const_global_mesh.py    # Global mesh parameters
│   ├── const_local_mesh.py     # Local mesh parameters
│   ├── material_property.py    # Material properties
│   └── simulation_params.py    # Simulation settings
└── utils/
    ├── logger.py           # Logging utilities
    └── step2str.py         # Step number formatting
```

## Module Descriptions

### 1. global_mesh.py - GlobalMesh Class
Generates global mesh with:
- Tensor product B-spline/NURBS mesh
- Graded refinement (m0 levels)
- Fine elements near crack, coarser elements far away
- Outputs: node.g.dat, elem.g.dat, weights.g.dat, index.g.dat

**Key Features:**
- Adjustable refinement levels (m0)
- Ratio between global and local element size (rGL)
- Symmetry considerations for efficient modeling

### 2. local_mesh.py - LocalMesh Class
Creates refined mesh around crack:
- Polar-like coordinate system
- Fine elements on crack surface (aL elements)
- Elements on ligament (lL elements)
- Through-thickness elements (HL elements)
- Outputs: node.l.dat, elem.l.dat, weights.l.dat

**Key Features:**
- Circular crack geometry
- Periodic boundary in angular direction
- Compatible with global mesh at interface

### 3. boundary.py - Boundary Class
Defines boundary conditions:
- Symmetry conditions on X=0, Y=-thi/2, Z=0 planes
- Remote stress application
- Interface constraints between global and local meshes
- Outputs: bc.g.dat, bc.l.dat, load.dat

**Key Features:**
- Automatic symmetry plane detection
- Distributed load calculation
- Separate BC files for global and local meshes

### 4. input_generator.py
Generates input.dat file containing:
- Problem type (DYNAMIC/STATIC)
- Mesh file references
- Material properties
- Time integration parameters (HHT method)
- Rayleigh damping coefficients
- Analysis settings
- Restart file references

### 5. initial.py
Handles restart for crack propagation:
- Reads previous step displacement, velocity, acceleration
- Initializes current step
- Manages dynamic crack growth
- Outputs: delta_u_init.dat, velocity_init.dat, acceleration_init.dat

### 6. linux_command.py
Executes SFEM solver:
- Runs sfem_linear executable
- Sets OpenMP threads
- Captures output
- Checks results
- Extracts and organizes output files

## Configuration Parameters

### const/const_global_mesh.py
```python
m0 = 4      # Number of mesh refinement levels
nx1 = 20    # Number of elements in X direction
ny1 = 2     # Number of elements in Y direction
rGL = 6     # Ratio between global and local element size
```

### const/const_local_mesh.py
```python
hL = 0.05e-3    # Local element size [m]
aL = 12         # Elements on crack surface
lL = 15         # Elements on ligament
HL = 8          # Elements through thickness
```

### const/material_property.py
```python
SigmaInfinity = 1.0e11  # Applied stress [Pa]
Nu = 0.3                # Poisson's ratio
EE = 2.06e11            # Young's modulus [Pa]
Rho = 7800.0            # Density [kg/m³]
```

### const/simulation_params.py
```python
c = 10.0e-3         # Crack radius [m]
thi = 1.0           # Thickness [m]
V = 400             # Crack velocity [m/s]
step_start = 0      # Starting step
step_end = 6        # Ending step
Alpha = 0.0         # HHT parameter
Beta = 0.25         # Newmark parameter
Gamma = 0.5         # Newmark parameter
OPENMP = 8          # OpenMP threads
```

## Usage

### Basic Execution
```python
python main.py
```

This will:
1. Generate meshes for each step (step_start to step_end)
2. Apply boundary conditions
3. Create input files
4. Execute SFEM solver
5. Handle restart between steps

### Step-by-step Execution
```python
from circular_crack.main import makemeshs
from circular_crack import linux_command

# Generate mesh for step 0
step = 0
REstart = 0  # 0 for first step, 1 for restart
l, g = makemeshs(step, REstart)

# Run solver
linux_command.run(step)
```

### Custom Parameters
Modify parameters in `const/` files before running:
```python
# Edit const/simulation_params.py
step_start = 0
step_end = 10
OPENMP = 16
```

## Output Files

### Mesh Files (per step)
- `node.g.dat` - Global mesh nodes (ID, X, Y, Z)
- `elem.g.dat` - Global mesh elements (connectivity)
- `weights.g.dat` - NURBS weights
- `node.l.dat` - Local mesh nodes
- `elem.l.dat` - Local mesh elements

### Boundary Condition Files
- `bc.g.dat` - Global mesh BCs (NodeID, DOF, Value)
- `bc.l.dat` - Local mesh BCs
- `load.dat` - Applied loads (NodeID, DOF, Force)

### Solver Input
- `input.dat` - Complete solver input file

### Results (generated by solver)
- `delta_u.dat` - Displacement field
- `velocity.dat` - Velocity field (dynamic)
- `acceleration.dat` - Acceleration field (dynamic)
- `stress.dat` - Stress field
- `strain.dat` - Strain field

## Workflow

1. **Step 0 (Initial)**
   - Generate initial meshes
   - Apply boundary conditions
   - Run solver with zero initial conditions

2. **Step 1 to N (Crack Propagation)**
   - Generate mesh with updated crack length
   - Load previous step results
   - Apply boundary conditions
   - Run solver with restart
   - Extract results

3. **Post-processing**
   - Results stored in `outputfolder/stepXXXXX/`
   - Can be visualized using Paraview or similar tools

## Notes

### Mesh Quality
- Element size ratio (rGL) affects solution accuracy
- Local mesh should be fine enough to capture crack tip fields
- Global mesh grading prevents excessive DOF

### Time Integration
- HHT-α method for dynamic analysis
- Rayleigh damping for numerical stability
- Time step determined by crack velocity

### Restart Strategy
- Each step uses previous displacement/velocity/acceleration
- Crack advances by c/stepall per step
- Meshes regenerated each step to follow crack growth

## Conversion from Mathematica

This Python implementation faithfully reproduces the Mathematica workflow:

1. **Mesh Generation**: Tensor product B-spline basis
2. **Geometry**: Circular crack with radial symmetry
3. **Refinement**: Multiple levels of graded mesh
4. **Boundary Conditions**: Symmetry planes and remote stress
5. **Time Integration**: HHT method with Rayleigh damping

Key differences:
- NumPy arrays instead of Mathematica lists
- Object-oriented design for better modularity
- Integrated logging for debugging
- Automatic directory management

## Troubleshooting

### "SFEM executable not found"
Compile the solver first:
```bash
cd ../sfem_linear
make
```

### "Missing output files"
Check solver log:
```bash
cat inputfiles/stepXXXXX/solver.log
```

### Mesh generation issues
Verify parameters in `const/` files:
- Element sizes must be positive
- Refinement levels reasonable (m0 < 10)
- Crack length < total domain size

## Future Enhancements

- [ ] Adaptive mesh refinement
- [ ] J-integral calculation
- [ ] Parallel step execution
- [ ] Visualization tools
- [ ] Result post-processing utilities
- [ ] Validation against analytical solutions

## References

Based on:
- Mathematica notebooks for circular crack analysis
- S-IGA method for crack propagation
- SFEM linear solver documentation

## License

Follow the license of the parent SFEM project.
