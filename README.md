# hS-IGA: 3D circular-crack implementation

This repository contains the three-dimensional circular-crack implementation used for the hS-IGA CMAME paper. Python code generates the B-spline-based global mesh, Lagrangian-type local mesh (near the crack front to ensure high accuracy), boundary conditions, and solver input. The accompanying sfem_linear Fortran solver performs the linear analysis. Static and dynamic modes, field output, and J-integral/DSIF post-processing are available.

The detailed method, validation, and numerical results are documented in the paper; this README only records how to install and run the released 3D code.

## Requirements

- Linux environment; the tested installation path is Ubuntu 22.04 or a compatible distribution.
- Python 3.10.
- Git, Make, a C/C++ and Fortran toolchain, CMake, and OpenMPI when building the solver.
- Authorised access to the collaborator-managed sfem_linear and nested Monolis repositories.

## Installation

The main repository records the exact sfem_linear commit required for this release. The solver is a private collaborator dependency, so first obtain access from the maintainers (or an approved mirror containing the recorded commit).

Clone the main repository normally, then run the setup helper:

    git clone <REPOSITORY-URL>
    cd <REPOSITORY-DIRECTORY>
    ./setup.sh --install-system-deps

If the listed system packages are already installed, omit --install-system-deps. The script creates a project-local .venv, installs requirements.txt, initialises sfem_linear and its nested submodules at the commit pinned by this repository, applies the project compatibility patch when needed, and builds the solver.

The script does not change system-wide compiler alternatives, build Python from source, edit shell startup files, or automatically advance sfem_linear to a later branch tip. To use an approved mirror or a local solver clone, set its URL explicitly:

    SFEM_LINEAR_REPO=<APPROVED-SOLVER-URL-OR-PATH> ./setup.sh

The supplied solver location must contain the sfem_linear commit recorded in this repository. Confirm the installation without modifying it with:

    ./setup.sh --check

## Representative commands

Run commands from circular_crack/ with the project Python environment.

    cd circular_crack

    # Inspect all supported controls.
    ../.venv/bin/python main.py --help

    # Generate input for one dynamic step without running the Fortran solver.
    ../.venv/bin/python main.py --meshonly --step_start 0 --step_end 1

    # Run one static validation analysis (requires a built solver).
    ../.venv/bin/python main.py --static_only

    # Post-process existing results.
    ../.venv/bin/python main.py --is_K --step_start 1 --step_end 10

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
