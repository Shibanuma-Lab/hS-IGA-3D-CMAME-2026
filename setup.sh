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
SOLVER_DIR=hs_iga
SOLVER_BIN="$SOLVER_DIR/bin/hs_iga"
MONOLIS_DIR="$SOLVER_DIR/submodule/monolis"
MONOLIS_LIB="$MONOLIS_DIR/lib/libmonolis_solver.a"
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
  --skip-solver          Do not initialise or build hs_iga.
  --force-rebuild        Rebuild Monolis and hs_iga.
  --check                Check an existing installation without modifying it.
  -h, --help             Show this help message.

Environment variables:
  PYTHON=/path/python3.10       Select the Python 3.10 interpreter.
  VENV_DIR=/path/to/.venv       Select a project virtual-environment directory.
  HS_IGA_REPO=<URL>        Use an alternate public hs_iga mirror or local clone.

The solver is pinned to the hs_iga commit recorded by this repository.
HS_IGA_REPO must provide that exact commit. It is not automatically
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
    commit=$(git rev-parse ":$SOLVER_DIR")
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
        libblas-dev
        liblapack-dev
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
    source_url=$(git config --file .gitmodules --get submodule.hs_iga.url || true)
    [[ -n "$source_url" ]] || die "No URL is configured for the hs_iga submodule."

    if [[ -v HS_IGA_REPO ]]; then
        source_url="$HS_IGA_REPO"
    fi

    if [[ -e "$SOLVER_DIR" ]] && ! is_solver_checkout; then
        # A normal clone can leave an empty directory for an uninitialised submodule.
        # Remove only that empty placeholder; preserve any non-empty user directory.
        if [[ -d "$SOLVER_DIR" ]] && rmdir "$SOLVER_DIR"; then
            info "Removed empty hs_iga submodule placeholder."
        else
            die "$SOLVER_DIR exists but is not a Git checkout. Move it aside manually, then rerun setup."
        fi
    fi

    if is_solver_checkout; then
        if [[ -n "$(git -C "$SOLVER_DIR" status --porcelain)" ]]; then
            die "$SOLVER_DIR has local changes. Commit, stash, or use a clean checkout before setup."
        fi
    fi

    info "Initialising hs_iga at recorded commit $expected_commit"
    git submodule sync -- "$SOLVER_DIR"

    # The override is local configuration only: it never rewrites .gitmodules.
    git config submodule.hs_iga.url "$source_url"

    if ! git submodule update --init --checkout "$SOLVER_DIR"; then
        die "Could not obtain hs_iga. Verify HTTPS access to the public repository or set HS_IGA_REPO to an alternate public mirror/local clone containing commit $expected_commit."
    fi

    actual_commit=$(git -C "$SOLVER_DIR" rev-parse HEAD)
    if [[ "$actual_commit" != "$expected_commit" ]]; then
        die "hs_iga commit mismatch. Expected $expected_commit but checked out $actual_commit."
    fi

    info "Initialising nested Monolis dependencies."
    if ! git -C "$SOLVER_DIR" submodule sync --recursive; then
        die "Could not synchronise nested hs_iga submodule URLs."
    fi
    # The public hS IGA manifest records Monolis with an SSH URL.  Use a local
    # HTTPS override so no GitHub account or SSH key is required.
    git -C "$SOLVER_DIR" config submodule.submodule/monolis.url "https://github.com/nqomorita/monolis.git"
    if ! git -C "$SOLVER_DIR" submodule update --init --recursive; then
        die "Could not initialise nested Monolis dependencies. Verify network access to their public Git remotes."
    fi

    [[ -d "$MONOLIS_DIR" ]] || die "Monolis was not initialised at $MONOLIS_DIR."
}

confirm_public_monolis_layout() {
    info "Using the Monolis revision pinned by public hs_iga; no private compatibility patch is required."
}

monolis_library_exists() {
    [[ -f "$MONOLIS_LIB" ]] || return 1
    return 0
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
    if ! monolis_library_exists; then
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
    require_command cmake

    [[ -x "$SOLVER_DIR/install_lib.sh" ]] || die "Missing public hs_iga dependency installer: $SOLVER_DIR/install_lib.sh"
    info "Building the Monolis dependencies pinned by public hs_iga."
    (
        cd "$SOLVER_DIR"
        ./install_lib.sh
    )

    [[ -f "$MONOLIS_LIB" ]] || die "Public hs_iga dependency build did not produce $MONOLIS_LIB."
}

build_solver() {
    require_command mpif90
    require_command mpicc

    info "Building hs_iga."
    make -C "$SOLVER_DIR" clean
    make -C "$SOLVER_DIR" FC="mpif90 -fopenmp" CC="mpicc -std=c99"

    [[ -x "$SOLVER_BIN" ]] || die "hs_iga build completed without producing $SOLVER_BIN."
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
        warn "hs_iga is not initialised."
        failures=1
    else
        actual_commit=$(git -C "$SOLVER_DIR" rev-parse HEAD)
        if [[ "$actual_commit" != "$expected_commit" ]]; then
            warn "hs_iga is at $actual_commit; expected $expected_commit."
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
    confirm_public_monolis_layout

    if monolis_needs_rebuild; then
        build_monolis
    else
        info "Monolis library is current."
    fi

    if solver_needs_rebuild; then
        build_solver
    else
        info "hs_iga binary is current."
    fi
else
    info "Skipping hs_iga initialisation and build."
fi

mkdir -p circular_crack/logs

info "Setup completed."
info "Run a representative command with:"
printf '  %s circular_crack/main.py --help\n' "$VENV_DIR/bin/python"
