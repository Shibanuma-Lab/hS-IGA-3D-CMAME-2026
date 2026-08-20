# hS-IGA: 3D circular-crack implementation

This repository contains the three-dimensional circular-crack implementation used for the hS-IGA CMAME paper. Python code generates the B-spline-based global mesh, Lagrangian-type local mesh (near the crack front to ensure high accuracy), boundary conditions, and solver input. The accompanying hs_iga Fortran solver performs the linear analysis. Static and dynamic modes, field output, and J-integral/DSIF post-processing are available.

The detailed method, validation, and numerical results are documented in the paper; this README only records how to install and run the released 3D code.

## Requirements

- Ubuntu 22.04/24.04 or Ubuntu under WSL 2. The setup script installs Python 3.10, the compiler toolchain, MPI, CMake, and BLAS/LAPACK development libraries when they are absent.
- Internet access to public GitHub and GitLab repositories, plus sudo permission on a new machine. No SSH key, account, or collaborator permission is required.

No prior Python, pip, virtual-environment, compiler, or MPI setup is required for the supported Ubuntu path.

## Installation

The main repository pins the exact public hs_iga commit required for this release.

On a brand-new Ubuntu or WSL installation, install Git once so that the repository can be obtained:

    sudo apt-get update
    sudo apt-get install -y git

Then clone the main repository and run the setup helper. No manual Python installation is required:

    git clone https://github.com/Shibanuma-Lab/hS-IGA-3D-CMAME-2026.git
    cd hS-IGA-3D-CMAME-2026
    ./setup.sh

The script detects missing Ubuntu dependencies and installs them automatically, including Python 3.10. If the configured Ubuntu repositories do not provide Python 3.10, it adds the deadsnakes PPA and continues. It then creates a project-local .venv, installs requirements.txt, initialises the pinned public hs_iga and Monolis submodules recursively over HTTPS, and builds the solver.

The public upstream manifest for Monolis uses an SSH URL; setup applies a local HTTPS override without modifying that upstream manifest, so no SSH key is needed. To prevent system package installation on an administrator-managed system, use ./setup.sh --skip-system-deps. To force a refresh of the supported Ubuntu packages, use ./setup.sh --install-system-deps.

The script does not change system-wide compiler alternatives, build Python from source, edit shell startup files, or automatically advance the pinned public solver commit. To use a public mirror or local public checkout, set its URL explicitly:

    HS_IGA_REPO=<PUBLIC-HS-IGA-URL-OR-PATH> ./setup.sh

The supplied public solver location must contain the hs_iga commit recorded in this repository. Confirm the installation without modifying it with:

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

    # Run the v=500 baseline dynamic case for steps 0 and 1 only.
    ../.venv/bin/python param_sweep_dynamic.py --velocities 500 --steps 2 --only-baseline --no-postprocess
    # Post-process existing results.
    ../.venv/bin/python main.py --is_K --step_start 1 --step_end 10

The normal Python launcher runs the public solver through one MPI process. For a dynamic analysis, omit --meshonly and select the step interval with --step_start and --step_end. The latter is exclusive; for example, --step_start 0 --step_end 10 processes steps 0 through 9. A restart assumes the previous-step files exist. Use --no_restart only when a fresh start is intended.

## Configuration and output

Edit the modules in circular_crack/const/ before a production run:

- simulation_params.py controls crack geometry, velocity, step range, HHT-alpha parameters, the default blended mass-lumping factor, and solver settings.
- material_property.py defines the elastic material parameters.
- const_global_mesh.py and const_local_mesh.py define the global and local discretisations.
- const_jintegral.py defines the J-integral post-processing domain.

The default dynamic launcher uses SFEM_MASS_LUMPING_ALPHA=0.02. Generated input is placed in circular_crack/inputfiles/. Solver results and logs are placed below circular_crack/results/ and circular_crack/logs/. These generated directories are intentionally ignored by Git.

## Scope and limitations

- The code models linear elastic, prescribed circular-crack configurations. The released examples do not autonomously predict initiation, crack-front direction, propagation velocity, or arrest.
- The bundled Sneddon data support the analytical static reference. Some optional comparison and sweep scripts additionally use case-specific FEM reference data under circular_crack/data/; confirm that the required data are present for the script you choose.
- Full simulations can be computationally demanding. Start with the one-step mesh-only command and record the configuration, compiler, MPI/OpenMP settings, and thread count for reproducibility.

## License

This code is released under the [MIT License](LICENSE). The public hs_iga solver and its public dependencies retain their own licenses.
