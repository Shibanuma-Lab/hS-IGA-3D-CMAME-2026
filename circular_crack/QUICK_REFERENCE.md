# Quick Reference Guide

## Installation & Setup

```bash
# Navigate to project directory
cd ~/S-IGA-circular-crack-in-3D-solid-linear/circular_crack

# Install Python dependencies
pip install numpy logzero

# Setup directories
python quickstart.py --setup

# Verify installation
python quickstart.py --config
```

## Running Simulations

### 1. Test Mode (No Solver)
```bash
# Run all tests
./test_mesh_generation.py

# Or using quickstart
./quickstart.py --test
```

### 2. Single Step
```bash
# Run step 0 only
./quickstart.py --step 0

# Or directly in Python
python -c "from main import makemeshs; from linux_command import run; makemeshs(0, 0); run(0)"
```

### 3. Full Simulation
```bash
# Run all steps (step_start to step_end)
./quickstart.py --run

# Or
python main.py
```

## Configuration

### Quick Edit Parameters
```python
# Edit const/simulation_params.py
step_start = 0      # Starting step
step_end = 10       # Ending step
OPENMP = 8          # Number of threads

# Edit const/const_global_mesh.py
rGL = 6             # Global/local ratio
m0 = 4              # Refinement levels

# Edit const/const_local_mesh.py
hL = 0.05e-3        # Element size [m]
aL = 12             # Crack elements
```

## File Locations

### Input Files
```
inputfiles/step00000/
  ├── node.g.dat      # Global nodes
  ├── elem.g.dat      # Global elements
  ├── bc.g.dat        # Boundary conditions
  ├── load.dat        # Applied loads
  └── input.dat       # Solver input
```

### Output Files
```
outputfolder/step00000/
  ├── delta_u.dat     # Displacement
  ├── stress.dat      # Stress field
  └── strain.dat      # Strain field
```

## Common Commands

```bash
# Show configuration
./quickstart.py --config

# Test mesh generation
./test_mesh_generation.py

# Run specific step
./quickstart.py --step 3

# Full simulation
./quickstart.py --run

# View logs
cat logs/arrest.log

# Check results
ls -lh outputfolder/step*/
```

## Python API Usage

```python
# Import modules
from global_mesh import GlobalMesh
from local_mesh import LocalMesh
from boundary import Boundary
import input_generator

# Generate global mesh
g = GlobalMesh(step=0)
g.make_global_mesh()
g.generate()

# Generate local mesh
l = LocalMesh(step=0)
l.make_local_mesh()
l.generate()

# Define boundary conditions
b = Boundary(l, g)
b.define_boundary(l, g)
b.generate()

# Generate input file
input_generator.generate(step=0, REstart=0)
```

## Parameter Ranges

| Parameter | Typical Range | Description |
|-----------|---------------|-------------|
| rGL       | 4-10          | Global/local element ratio |
| m0        | 3-5           | Refinement levels |
| hL        | 0.01-0.1 mm   | Local element size |
| aL        | 10-20         | Crack surface elements |
| lL        | 10-20         | Ligament elements |
| HL        | 5-15          | Through-thickness elements |
| OPENMP    | 4-16          | OpenMP threads |

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Solver not found | `cd ../sfem_linear && make` |
| Import errors | `pip install numpy logzero` |
| Missing directories | `./quickstart.py --setup` |
| Permission denied | `chmod +x *.py` |
| Memory error | Reduce mesh density (increase rGL) |

## Performance Tips

1. **Optimal Threading**: Set OPENMP to physical core count
2. **Mesh Efficiency**: Use rGL=6-8 for balanced accuracy/speed
3. **Storage**: Clean old results to save disk space
4. **Monitoring**: Check logs/arrest.log for progress

## Example Workflows

### Quick Test Run
```bash
# 1. Setup
./quickstart.py --setup

# 2. Test mesh generation
./quickstart.py --test

# 3. Run step 0
./quickstart.py --step 0

# 4. Check results
ls -lh inputfiles/step00000/
```

### Production Run
```bash
# 1. Configure parameters
nano const/simulation_params.py  # Edit step_end=20, OPENMP=16

# 2. Run simulation
./quickstart.py --run

# 3. Monitor progress
tail -f logs/arrest.log

# 4. Extract results
ls -lh outputfolder/step*/
```

### Parameter Study
```python
# Create script: param_study.py
import sys
from const import const_global_mesh as cgm
from main import makemeshs
from linux_command import run

for rGL in [4, 6, 8, 10]:
    cgm.rGL = rGL
    print(f"Running with rGL={rGL}")
    l, g = makemeshs(step=0, REstart=0)
    run(0)
```

## File Format Reference

### node.g.dat
```
NodeID  X-coord         Y-coord         Z-coord
1       0.000000e+00    -5.000000e-01   0.000000e+00
2       5.000000e-04    -5.000000e-01   0.000000e+00
```

### elem.g.dat
```
ElemID  Node1  Node2  ...  Node64
1       1      2      ...  64
2       2      3      ...  65
```

### bc.g.dat
```
NodeID  DOF  Value
1       1    0.0
2       2    0.0
```

### load.dat
```
NodeID  DOF  Force
100     1    1.5e6
101     1    1.5e6
```

## Key Classes & Functions

| Module | Class/Function | Purpose |
|--------|----------------|---------|
| global_mesh | GlobalMesh | Generate global mesh |
| local_mesh | LocalMesh | Generate local mesh |
| boundary | Boundary | Define BCs and loads |
| input_generator | generate() | Create input.dat |
| initial | initial() | Handle restart |
| linux_command | run() | Execute solver |

## Important Notes

1. **File Naming**: Uses 5-digit step numbers (00000, 00001, etc.)
2. **Indexing**: 1-based indexing for compatibility with Fortran
3. **Units**: SI units (meters, Pascals, kg)
4. **Restart**: Automatically handles step > 0
5. **Symmetry**: Uses quarter or eighth model with symmetry BCs

## Getting Help

1. **Documentation**: See README.md and IMPLEMENTATION_SUMMARY.md
2. **Code Comments**: All functions have docstrings
3. **Test Suite**: test_mesh_generation.py shows usage examples
4. **Workflow**: See WORKFLOW_DIAGRAM.py for visual guide

## Contact Information

For questions or issues:
- Check code documentation
- Review test examples
- Refer to Mathematica notebooks for algorithm details
