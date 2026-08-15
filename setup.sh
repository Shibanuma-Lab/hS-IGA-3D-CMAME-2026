#!/usr/bin/env bash
# ============================================================================
# hS-IGA 3D circular-crack implementation: local installation helper
#
# This script intentionally keeps all project dependencies inside the checkout.
# It does not change system-wide compiler alternatives, build Python from source,
# or edit the user's shell startup files.
# ============================================================================

set -Eeuo pipefail
IFS=$'\n\t'

SCRIPT_NAME=$(basename "$0")
PROJECT_ROOT=$(cd -- "$(dirname -- "$0")" && pwd)
SOLVER_DIR=sfem_linear
SOLVER_BIN="$SOLVER_DIR/bin/sfem_linear"
MONOLIS_DIR="$SOLVER_DIR/submodule/monolis"
MONOLIS_LIB="$MONOLIS_DIR/lib/libmonolis_solver.a"
MONOLIS_PATCH="$PROJECT_ROOT/patches/monolis-siga-atomic-openmp.patch"
if [[ -v VENV_DIR ]]; then
    VENV_DIR="$VENV_DIR"
else
    VENV_DIR="$PROJECT_ROOT/.venv"
fi
INSTALL_SYSTEM_DEPS=0
SKIP_PYTHON=0
SKIP_SOLVER=0
CHECK_ONLY=0
FORCE_REBUILD=0
PYTHON_CMD=

info() {
    printf '[INFO] %s\n' "$*"
}

warn() {
    printf '[WARN] %s\n' "$*" >&2
}

die() {
    printf '[ERROR] %s\n' "$*" >&2
    exit 1
}

usage() {
    cat <<'USAGE'
Usage: ./setup.sh [options]

Options:
  --install-system-deps  Install the required Ubuntu packages with apt-get.
  --skip-python          Do not create or update the project-local .venv.
  --skip-solver          Do not initialise or build sfem_linear.
  --force-rebuild        Rebuild Monolis and sfem_linear.
  --check                Check an existing installation without modifying it.
  -h, --help             Show this help message.

Environment variables:
  PYTHON=/path/python3.10       Select the Python 3.10 interpreter.
  VENV_DIR=/path/to/.venv       Select a project virtual-environment directory.
  SFEM_LINEAR_REPO=<URL>        Use an approved sfem_linear mirror or local clone.

The solver is pinned to the sfem_linear commit recorded by this repository.
SFEM_LINEAR_REPO must provide that exact commit. It is not automatically
fast-forwarded to a newer branch tip.
USAGE
}

while (( $# > 0 )); do
    case "$1" in
        --install-system-deps)
            INSTALL_SYSTEM_DEPS=1
            ;;
        --skip-python)
            SKIP_PYTHON=1
            ;;
        --skip-solver)
            SKIP_SOLVER=1
            ;;
        --force-rebuild)
            FORCE_REBUILD=1
            ;;
        --check)
            CHECK_ONLY=1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1 (run ./$SCRIPT_NAME --help)"
            ;;
    esac
    shift
done



cd "$PROJECT_ROOT"

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    die "Run this script from a Git checkout of the 3D repository."
fi

if [[ ! -d circular_crack ]]; then
    die "circular_crack/ is missing; this does not look like the 3D project root."
fi

if [[ ! -f requirements.txt ]]; then
    die "requirements.txt is missing."
fi

if [[ ! -f .gitmodules ]]; then
    die ".gitmodules is missing. Please update the main repository checkout before running setup."
fi

expected_solver_commit() {
    local commit
    commit=$(git ls-tree HEAD -- "$SOLVER_DIR" | awk '$1 == "160000" {print $3}')
    if [[ -z "$commit" ]]; then
        die "The main repository does not record a commit for $SOLVER_DIR."
    fi
    printf '%s\n' "$commit"
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

resolve_python() {
    local version

    if [[ -v PYTHON ]]; then
        PYTHON_CMD="$PYTHON"
    elif command -v python3.10 >/dev/null 2>&1; then
        PYTHON_CMD=$(command -v python3.10)
    else
        die "Python 3.10 is required. Install it first, or pass PYTHON=/path/to/python3.10."
    fi

    require_command "$PYTHON_CMD"
    version=$("$PYTHON_CMD" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
    if [[ "$version" != "3.10" ]]; then
        die "Python 3.10 is required by Pipfile.lock; selected interpreter is Python $version."
    fi
}

run_as_root() {
    if (( EUID == 0 )); then
        "$@"
    elif command -v sudo >/dev/null 2>&1; then
        sudo "$@"
    else
        die "Administrator privileges are required for --install-system-deps, but sudo is unavailable."
    fi
}

install_system_dependencies() {
    local packages

    require_command apt-get
    packages=(
        build-essential
        cmake
        git
        make
        pkg-config
        gfortran
        openmpi-bin
        libopenmpi-dev
        python3.10
        python3.10-venv
        ca-certificates
    )

    info "Installing required Ubuntu packages (no system upgrade and no compiler-alternative changes)."
    run_as_root apt-get update
    run_as_root env DEBIAN_FRONTEND=noninteractive apt-get install -y "$packages"
}

create_python_environment() {
    local venv_version

    resolve_python

    if [[ -e "$VENV_DIR" && ! -x "$VENV_DIR/bin/python" ]]; then
        die "$VENV_DIR exists but is not a usable virtual environment. Choose another VENV_DIR or remove it manually."
    fi

    if [[ ! -x "$VENV_DIR/bin/python" ]]; then
        info "Creating project-local Python environment: $VENV_DIR"
        "$PYTHON_CMD" -m venv "$VENV_DIR"
    fi

    venv_version=$("$VENV_DIR/bin/python" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
    if [[ "$venv_version" != "3.10" ]]; then
        die "$VENV_DIR uses Python $venv_version, but this release requires Python 3.10. Choose another VENV_DIR or recreate it manually."
    fi

    info "Installing Python packages into $VENV_DIR"
    "$VENV_DIR/bin/python" -m pip install --upgrade pip
    "$VENV_DIR/bin/python" -m pip install -r requirements.txt

    "$VENV_DIR/bin/python" - <<'PY'
import importlib

modules = ("numpy", "scipy", "logzero", "pandas", "openpyxl", "numba", "matplotlib")
missing = []
for name in modules:
    try:
        importlib.import_module(name)
    except Exception as exc:
        missing.append(f"{name}: {exc}")

if missing:
    raise SystemExit("Invalid Python environment:\n" + "\n".join(missing))
PY

    info "Python environment is ready."
}

is_solver_checkout() {
    [[ -e "$SOLVER_DIR/.git" ]] && git -C "$SOLVER_DIR" rev-parse --is-inside-work-tree >/dev/null 2>&1
}

initialise_solver_checkout() {
    local expected_commit
    local source_url
    local actual_commit

    expected_commit=$(expected_solver_commit)
    source_url=$(git config --file .gitmodules --get submodule.sfem_linear.url || true)
    [[ -n "$source_url" ]] || die "No URL is configured for the sfem_linear submodule."

    if [[ -v SFEM_LINEAR_REPO ]]; then
        source_url="$SFEM_LINEAR_REPO"
    fi

    if [[ -e "$SOLVER_DIR" ]] && ! is_solver_checkout; then
        # A normal clone can leave an empty directory for an uninitialised submodule.
        # Remove only that empty placeholder; preserve any non-empty user directory.
        if [[ -d "$SOLVER_DIR" ]] && rmdir "$SOLVER_DIR"; then
            info "Removed empty sfem_linear submodule placeholder."
        else
            die "$SOLVER_DIR exists but is not a Git checkout. Move it aside manually, then rerun setup."
        fi
    fi

    if is_solver_checkout; then
        if [[ -n "$(git -C "$SOLVER_DIR" status --porcelain)" ]]; then
            die "$SOLVER_DIR has local changes. Commit, stash, or use a clean checkout before setup."
        fi
    fi

    info "Initialising sfem_linear at recorded commit $expected_commit"
    git submodule sync -- "$SOLVER_DIR"

    # The override is local configuration only: it never rewrites .gitmodules.
    git config submodule.sfem_linear.url "$source_url"

    if ! git submodule update --init --checkout "$SOLVER_DIR"; then
        die "Could not obtain sfem_linear. This solver is collaborator-managed; request access or set SFEM_LINEAR_REPO to an approved mirror/local clone containing commit $expected_commit."
    fi

    actual_commit=$(git -C "$SOLVER_DIR" rev-parse HEAD)
    if [[ "$actual_commit" != "$expected_commit" ]]; then
        die "sfem_linear commit mismatch. Expected $expected_commit but checked out $actual_commit."
    fi

    info "Initialising nested Monolis dependencies."
    if ! git -C "$SOLVER_DIR" submodule sync --recursive; then
        die "Could not synchronise nested sfem_linear submodule URLs."
    fi
    if ! git -C "$SOLVER_DIR" submodule update --init --recursive; then
        die "Could not initialise nested Monolis dependencies. Verify access to their Git remotes."
    fi

    [[ -d "$MONOLIS_DIR" ]] || die "Monolis was not initialised at $MONOLIS_DIR."
}

monolis_atomic_source() {
    local candidate
    for candidate in "$MONOLIS_DIR/src/matrix/spmat_handler.f90" "$MONOLIS_DIR/src/matrix/sparse_util.f90"; do
        if [[ -f "$candidate" ]] && grep -q "subroutine monolis_add_scalar_to_sparse_matrix_atomic" "$candidate"; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

monolis_is_compatible() {
    monolis_atomic_source >/dev/null 2>&1 && grep -q -- "-fopenmp" "$MONOLIS_DIR/Makefile"
}

ensure_monolis_compatibility() {
    if monolis_is_compatible; then
        info "Monolis already has the required atomic assembly and OpenMP support."
        return
    fi

    [[ -f "$MONOLIS_PATCH" ]] || die "Missing project compatibility patch: $MONOLIS_PATCH"

    if git -C "$MONOLIS_DIR" apply --check "$MONOLIS_PATCH"; then
        info "Applying the project-local Monolis compatibility patch."
        git -C "$MONOLIS_DIR" apply "$MONOLIS_PATCH"
    elif git -C "$MONOLIS_DIR" apply --reverse --check "$MONOLIS_PATCH"; then
        info "The project-local Monolis compatibility patch is already applied."
    else
        die "The Monolis source is incompatible with the project patch. Update the pinned solver dependency rather than forcing this installation."
    fi

    monolis_is_compatible || die "Monolis compatibility checks failed after applying the project patch."
}

library_has_atomic_symbol() {
    [[ -f "$MONOLIS_LIB" ]] || return 1
    nm -a "$MONOLIS_LIB" 2>/dev/null | grep -F "monolis_add_scalar_to_sparse_matrix_atomic" >/dev/null
}
monolis_dependencies_ready() {
    local dependency

    for dependency in lib/libmetis.a lib/libparmetis.a lib/libmonolis_utils.a lib/libgedatsu.a lib/libggtools.a; do
        [[ -f "$MONOLIS_DIR/$dependency" ]] || return 1
    done
}

monolis_needs_rebuild() {
    if (( FORCE_REBUILD == 1 )); then
        return 0
    fi
    if ! monolis_dependencies_ready; then
        return 0
    fi
    if [[ ! -f "$MONOLIS_LIB" ]]; then
        return 0
    fi
    if ! library_has_atomic_symbol; then
        return 0
    fi
    find "$MONOLIS_DIR/src" "$MONOLIS_DIR/Makefile" -type f -newer "$MONOLIS_LIB" -print -quit | grep -q .
}

solver_needs_rebuild() {
    if (( FORCE_REBUILD == 1 )); then
        return 0
    fi
    if [[ ! -x "$SOLVER_BIN" ]]; then
        return 0
    fi
    if [[ "$MONOLIS_LIB" -nt "$SOLVER_BIN" ]]; then
        return 0
    fi
    find "$SOLVER_DIR/src" "$SOLVER_DIR/Makefile" -type f -newer "$SOLVER_BIN" -print -quit | grep -q .
}

build_monolis() {
    require_command mpif90
    require_command mpicc
    require_command nm

    info "Building Monolis with MPI, METIS, and OpenMP support."
    if ! monolis_dependencies_ready; then
        [[ -x "$MONOLIS_DIR/install_lib.sh" ]] || die "Missing Monolis dependency installer: $MONOLIS_DIR/install_lib.sh"
        info "Building the nested Monolis libraries required by sfem_linear."
        (
            cd "$MONOLIS_DIR"
            ./install_lib.sh METIS
        )
    fi
    make -C "$MONOLIS_DIR" clean >/dev/null 2>&1 || true
    make -C "$MONOLIS_DIR" FLAGS=MPI,METIS FC=mpif90 CC=mpicc

    library_has_atomic_symbol || die "The Monolis library was built without the required atomic sparse-matrix symbol."
}

build_solver() {
    require_command mpif90
    require_command mpicc

    info "Building sfem_linear."
    make -C "$SOLVER_DIR" clean
    make -C "$SOLVER_DIR" FC="mpif90 -fopenmp" CC="mpicc -std=c99"

    [[ -x "$SOLVER_BIN" ]] || die "sfem_linear build completed without producing $SOLVER_BIN."
}

check_installation() {
    local expected_commit
    local actual_commit
    local failures=0

    expected_commit=$(expected_solver_commit)

    if [[ ! -x "$VENV_DIR/bin/python" ]]; then
        warn "Missing Python environment: $VENV_DIR"
        failures=1
    fi

    if ! is_solver_checkout; then
        warn "sfem_linear is not initialised."
        failures=1
    else
        actual_commit=$(git -C "$SOLVER_DIR" rev-parse HEAD)
        if [[ "$actual_commit" != "$expected_commit" ]]; then
            warn "sfem_linear is at $actual_commit; expected $expected_commit."
            failures=1
        fi
    fi

    if [[ ! -x "$SOLVER_BIN" ]]; then
        warn "Missing solver binary: $SOLVER_BIN"
        failures=1
    fi

    if (( failures != 0 )); then
        die "Installation check failed. Run ./$SCRIPT_NAME to install missing components."
    fi

    info "Installation check passed."
}

if (( CHECK_ONLY == 1 )); then
    check_installation
    exit 0
fi

if (( INSTALL_SYSTEM_DEPS == 1 )); then
    install_system_dependencies
fi

require_command git
require_command make

if (( SKIP_PYTHON == 0 )); then
    create_python_environment
else
    info "Skipping Python environment setup."
fi

if (( SKIP_SOLVER == 0 )); then
    initialise_solver_checkout
    ensure_monolis_compatibility

    if monolis_needs_rebuild; then
        build_monolis
    else
        info "Monolis library is current."
    fi

    if solver_needs_rebuild; then
        build_solver
    else
        info "sfem_linear binary is current."
    fi
else
    info "Skipping sfem_linear initialisation and build."
fi

mkdir -p circular_crack/logs

info "Setup completed."
info "Run a representative command with:"
printf '  %s circular_crack/main.py --help\n' "$VENV_DIR/bin/python"
