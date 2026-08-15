# hS-IGA: 3D circular-crack implementation

This repository contains the three-dimensional circular-crack implementation used for the hS-IGA CMAME paper. Python code generates the NURBS-based global mesh, local crack-front discretisation, boundary conditions, and solver input. The accompanying sfem_linear Fortran solver performs the linear analysis. Static and dynamic modes, field output, and J-integral/DSIF post-processing are available.

The detailed method, validation, and numerical results are documented in the paper; this README only records how to install and run the released 3D code.

## Requirements

- Linux environment (the solver build and execution scripts target Linux)
- Python 3.10 with the packages in requirements.txt
- A C/C++ and Fortran toolchain, CMake, Make, and OpenMPI when rebuilding the solver
- Git, including submodule support

The source includes Pipfile.lock for a reproducible Python environment.

## Installation

Clone the repository with its solver submodule, then run the setup script:

    git clone --recurse-submodules <PUBLIC-REPOSITORY-URL>
    cd <REPOSITORY-DIRECTORY>
    ./setup.sh

setup.sh installs the Python packages, initialises the nested solver dependencies, and builds sfem_linear when no compatible executable is present. It may request administrator privileges to install compiler and MPI packages. To use an existing Python environment instead, install the packages with:

    python -m pip install -r requirements.txt

Confirm that sfem_linear/bin/sfem_linear exists before attempting an analysis.

**Release requirement:** the current sfem_linear submodule URL is private. Before making this repository public, replace it with the approved public solver URL (and ensure its nested dependencies are accessible), or include the solver source in the public release. A public clone must be tested with git clone --recurse-submodules from an account without laboratory access.

## Representative commands

Run commands from circular_crack/.

    cd circular_crack

    # Inspect all supported controls.
    python main.py --help

    # Generate input for one dynamic step without running the Fortran solver.
    python main.py --meshonly --step_start 0 --step_end 1

    # Run one static validation analysis (requires a built solver).
    python main.py --static_only

    # Post-process existing results.
    python main.py --is_K --step_start 1 --step_end 10

For a normal dynamic analysis, omit --meshonly and select the step interval with --step_start and --step_end. The latter is exclusive; for example, --step_start 0 --step_end 10 processes steps 0 through 9. A restart assumes the previous-step files exist. Use --no_restart only when a fresh start is intended.

## Configuration and output

Edit the modules in circular_crack/const/ before a production run:

- simulation_params.py controls crack geometry, velocity, step range, and solver settings.
- material_property.py defines the elastic material parameters.
- const_global_mesh.py and const_local_mesh.py define the global and local discretisations.
- const_jintegral.py defines the J-integral post-processing domain.

Generated input is placed in circular_crack/inputfiles/. Solver results and logs are placed below circular_crack/results/ and circular_crack/logs/. These generated directories are intentionally ignored by Git.

## Scope and limitations

- The code models linear elastic, prescribed circular-crack configurations. The released examples do not autonomously predict initiation, crack-front direction, propagation velocity, or arrest.
- The bundled Sneddon data support the analytical static reference. Some optional comparison and sweep scripts additionally use case-specific FEM reference data under circular_crack/data/; confirm that the required data are present for the script you choose.
- Full simulations can be computationally demanding. Start with the one-step mesh-only command and record the configuration, compiler, and thread count for reproducibility.

## License

This code is released under the [MIT License](LICENSE).
