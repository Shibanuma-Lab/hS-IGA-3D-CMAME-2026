#!/bin/bash

# ============================================================================
# S-IGA Circular Crack 3D Solver - Installation Script
# ============================================================================
# This script automates the installation process for new computers
#
# Prerequisites:
#   - Git repository already cloned
#   - Run from project root directory
#
# Usage:
#   chmod +x setup.sh
#   ./setup.sh
#
# Optional environment variables:
#   SETUP_APT_UPGRADE=1       Run apt-get upgrade before installing packages
#   FORCE_REBUILD_SOLVER=1    Rebuild monolis and sfem_linear even if binary exists
#   SFEM_LINEAR_REPO=<url>    Override the sfem_linear clone URL
#   SFEM_LINEAR_BRANCH=<name> Override the sfem_linear branch (default: tianyu_IGA)
# ============================================================================

set -euo pipefail  # Exit on error, unset variables, and failed pipelines

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print functions
print_header() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

SUDO=""
if [ "$(id -u)" -ne 0 ]; then
    SUDO="sudo"
fi

APT_UPDATED=0

apt_update_once() {
    if [ "$APT_UPDATED" -eq 0 ]; then
        print_info "Running apt update..."
        $SUDO apt-get update
        APT_UPDATED=1
    fi
}

install_apt_packages() {
    if ! command -v apt-get &> /dev/null; then
        print_warning "apt-get not found. Skipping automatic system package installation."
        print_info "Please install build tools, GCC/G++/GFortran 11, OpenMPI, CMake, Git, and Wget manually."
        return
    fi

    if [ "${SETUP_APT_UPGRADE:-0}" = "1" ]; then
        apt_update_once
        print_info "Running apt upgrade because SETUP_APT_UPGRADE=1..."
        $SUDO apt-get upgrade -y
    fi

    local packages=(
        build-essential
        cmake
        make
        git
        wget
        ca-certificates
        tar
        pkg-config
        binutils
        gcc-11
        g++-11
        gfortran-11
        gfortran
        openmpi-doc
        openmpi-bin
        libopenmpi-dev
        zlib1g-dev
        libncurses5-dev
        libgdbm-dev
        libnss3-dev
        libssl-dev
        libsqlite3-dev
        libreadline-dev
        libffi-dev
        libbz2-dev
        liblzma-dev
        tk-dev
    )

    local missing=()
    local pkg
    for pkg in "${packages[@]}"; do
        if ! dpkg -s "$pkg" &> /dev/null; then
            missing+=("$pkg")
        fi
    done

    if [ "${#missing[@]}" -eq 0 ]; then
        print_success "System build dependencies are already installed"
        return
    fi

    apt_update_once
    print_info "Installing missing system packages: ${missing[*]}"
    $SUDO apt-get install -y "${missing[@]}"
}

register_compiler_alternatives() {
    local name="$1"
    local path
    local version
    local priority

    shopt -s nullglob
    for path in /usr/bin/"$name"-[0-9]*; do
        version="${path##*-}"
        if [[ "$version" =~ ^[0-9]+$ ]]; then
            priority=$((version * 10))
            $SUDO update-alternatives --install "/usr/bin/$name" "$name" "$path" "$priority"
        fi
    done
    shopt -u nullglob
}

configure_gcc_11() {
    local required=("/usr/bin/gcc-11" "/usr/bin/g++-11" "/usr/bin/gfortran-11")
    local exe

    for exe in "${required[@]}"; do
        if [ ! -x "$exe" ]; then
            print_error "$exe not found after package installation"
            exit 1
        fi
    done

    print_info "Registering compiler alternatives and selecting GCC/G++/GFortran 11..."
    register_compiler_alternatives "gcc"
    register_compiler_alternatives "g++"
    register_compiler_alternatives "gfortran"

    $SUDO update-alternatives --set gcc /usr/bin/gcc-11
    $SUDO update-alternatives --set g++ /usr/bin/g++-11
    $SUDO update-alternatives --set gfortran /usr/bin/gfortran-11

    export CC=/usr/bin/gcc-11
    export CXX=/usr/bin/g++-11
    export FC=/usr/bin/gfortran-11
    export OMPI_CC=/usr/bin/gcc-11
    export OMPI_CXX=/usr/bin/g++-11
    export OMPI_FC=/usr/bin/gfortran-11

    print_success "Compiler alternatives selected:"
    print_info "gcc: $(gcc -dumpfullversion -dumpversion)"
    print_info "g++: $(g++ -dumpfullversion -dumpversion)"
    print_info "gfortran: $(gfortran -dumpfullversion -dumpversion)"
}

ensure_system_dependencies() {
    print_header "Step 0: Installing System Dependencies"
    install_apt_packages
    configure_gcc_11
    echo ""
}

ensure_system_dependencies

# ============================================================================
# Step 1: Check directory structure
# ============================================================================
print_header "Step 1: Checking Directory Structure"

if [ ! -d "circular_crack" ]; then
    print_error "circular_crack directory not found!"
    print_info "Please ensure you are in the project root directory"
    exit 1
fi

print_success "circular_crack directory verified"

# Check if sfem_linear exists and is empty
if [ -d "sfem_linear" ]; then
    if [ -z "$(ls -A sfem_linear)" ]; then
        print_warning "sfem_linear directory is empty. Removing..."
        rmdir sfem_linear
        print_success "Empty sfem_linear directory removed"
    else
        print_info "sfem_linear directory exists with content"
    fi
fi

echo ""

# ============================================================================
# Step 2: Check/Install Python 3.10
# ============================================================================
print_header "Step 2: Checking Python 3.10"

PYTHON_VERSION="3.10.6"
PYTHON_EXEC="/usr/local/bin/python3.10"

if command -v python3.10 &> /dev/null; then
    CURRENT_VERSION=$(python3.10 --version | awk '{print $2}')
    print_success "Python 3.10 is already installed (version: $CURRENT_VERSION)"
    PYTHON_EXEC=$(command -v python3.10)
else
    print_warning "Python 3.10 not found. Installing..."

    TMP_BUILD_DIR=$(mktemp -d /tmp/siga-python-build.XXXXXX)
    print_info "Downloading Python $PYTHON_VERSION..."
    pushd "$TMP_BUILD_DIR" > /dev/null
    wget "https://www.python.org/ftp/python/$PYTHON_VERSION/Python-$PYTHON_VERSION.tgz"

    print_info "Extracting and building Python..."
    tar -xf "Python-$PYTHON_VERSION.tgz"
    pushd "Python-$PYTHON_VERSION" > /dev/null
    ./configure --enable-optimizations --with-ensurepip=install
    make -j"$(nproc)"

    print_info "Installing Python 3.10..."
    $SUDO make altinstall
    popd > /dev/null
    popd > /dev/null
    rm -rf "$TMP_BUILD_DIR"

    # Verify installation
    if python3.10 --version &> /dev/null; then
        print_success "Python 3.10 installed successfully!"
        PYTHON_EXEC=$(command -v python3.10)
    else
        print_error "Python 3.10 installation failed!"
        exit 1
    fi
fi

echo ""

# ============================================================================
# Step 3: Check/Install Pipenv
# ============================================================================
print_header "Step 3: Checking Pipenv"

if command -v pipenv &> /dev/null; then
    PIPENV_VERSION=$(pipenv --version | awk '{print $3}')
    print_success "Pipenv is already installed (version: $PIPENV_VERSION)"
else
    print_warning "Pipenv not found. Installing..."

    # Check if pip is available for Python 3.10
    if ! "$PYTHON_EXEC" -m pip --version &> /dev/null; then
        print_warning "pip not found for Python 3.10. Installing..."

        if "$PYTHON_EXEC" -m ensurepip --upgrade &> /dev/null; then
            print_success "pip installed via ensurepip"
        else
            print_info "Installing pip using get-pip.py..."
            TMP_PIP_DIR=$(mktemp -d /tmp/siga-pip-build.XXXXXX)
            pushd "$TMP_PIP_DIR" > /dev/null
            wget -q https://bootstrap.pypa.io/get-pip.py
            "$PYTHON_EXEC" get-pip.py --user
            popd > /dev/null
            rm -rf "$TMP_PIP_DIR"
        fi

        export PATH="$HOME/.local/bin:$PATH"

        if "$PYTHON_EXEC" -m pip --version &> /dev/null; then
            print_success "pip installed successfully!"
        else
            print_error "Failed to install pip for Python 3.10!"
            print_info "Please install pip manually, then re-run this script."
            exit 1
        fi
    fi

    print_info "Installing Pipenv..."
    "$PYTHON_EXEC" -m pip install --user pipenv

    # Add to PATH in .bashrc if not already present
    BASHRC="$HOME/.bashrc"
    PATH_EXPORT='export PATH="$HOME/.local/bin:$PATH"'

    if ! grep -q "$PATH_EXPORT" "$BASHRC" 2>/dev/null; then
        print_info "Adding Pipenv to PATH in .bashrc..."
        echo "" >> "$BASHRC"
        echo "# Added by S-IGA setup script" >> "$BASHRC"
        echo "$PATH_EXPORT" >> "$BASHRC"
        print_success "PATH updated in .bashrc"
    fi

    # Update current session
    export PATH="$HOME/.local/bin:$PATH"

    # Verify installation
    if command -v pipenv &> /dev/null; then
        print_success "Pipenv installed successfully!"
        print_warning "Note: You may need to restart your terminal or run: source ~/.bashrc"
    else
        print_error "Pipenv installation failed!"
        print_info "Please run: source ~/.bashrc"
        print_info "Then verify with: pipenv --version"
        exit 1
    fi
fi

echo ""

# ============================================================================
# Step 4: Setup Python Virtual Environment
# ============================================================================
print_header "Step 4: Setting Up Virtual Environment"

print_info "Activating Pipenv with Python 3.10..."
export PIPENV_VENV_IN_PROJECT=1  # Create .venv in project directory
pipenv --python "$PYTHON_EXEC"

print_info "Installing Python dependencies from Pipfile..."
if pipenv install; then
    print_success "Virtual environment created and dependencies installed!"
else
    print_error "Failed to install dependencies!"
    exit 1
fi

echo ""

# ============================================================================
# Step 5: Check/Create logs directory
# ============================================================================
print_header "Step 5: Checking Project Structure"

LOGS_DIR="circular_crack/logs"
if [ ! -d "$LOGS_DIR" ]; then
    print_warning "logs directory not found. Creating..."
    mkdir -p "$LOGS_DIR"
    print_success "Created $LOGS_DIR"
else
    print_success "logs directory exists"
fi

echo ""

# ============================================================================
# Step 6: Install and Build Fortran Solver
# ============================================================================
print_header "Step 6: Installing Fortran Solver"

SOLVER_DIR="sfem_linear"
SOLVER_BIN="$SOLVER_DIR/bin/sfem_linear"
SFEM_LINEAR_REPO="${SFEM_LINEAR_REPO:-git@gitlab.com:morita/sfem_linear.git}"
SFEM_LINEAR_BRANCH="${SFEM_LINEAR_BRANCH:-tianyu_IGA}"
MONOLIS_PATCH="$PROJECT_ROOT/patches/monolis-siga-atomic-openmp.patch"

prepare_sfem_linear_repo() {
    if [ ! -d "$SOLVER_DIR" ]; then
        print_warning "sfem_linear directory not found. Cloning from GitLab..."
        print_info "Cloning $SFEM_LINEAR_REPO..."
        if git clone "$SFEM_LINEAR_REPO" "$SOLVER_DIR"; then
            print_success "sfem_linear repository cloned successfully"
        else
            print_error "Failed to clone sfem_linear repository!"
            print_info "Set SFEM_LINEAR_REPO to another clone URL if this machine does not have GitLab SSH access."
            exit 1
        fi
    fi

    if git -C "$SOLVER_DIR" rev-parse --is-inside-work-tree &> /dev/null; then
        CURRENT_BRANCH=$(git -C "$SOLVER_DIR" branch --show-current 2>/dev/null || true)
        if [ "$CURRENT_BRANCH" != "$SFEM_LINEAR_BRANCH" ]; then
            print_info "Switching sfem_linear to branch: $SFEM_LINEAR_BRANCH"
            if git -C "$SOLVER_DIR" checkout "$SFEM_LINEAR_BRANCH"; then
                print_success "Switched to branch $SFEM_LINEAR_BRANCH"
            else
                print_error "Failed to checkout branch $SFEM_LINEAR_BRANCH"
                exit 1
            fi
        fi
    else
        print_error "$SOLVER_DIR exists but is not a git repository"
        exit 1
    fi
}

ensure_monolis_source() {
    local monolis_dir="$SOLVER_DIR/submodule/monolis"

    if [ ! -d "$monolis_dir" ] || [ -z "$(ls -A "$monolis_dir" 2>/dev/null)" ]; then
        print_info "Initializing sfem_linear submodules..."
        git -C "$SOLVER_DIR" submodule update --init --recursive
    fi

    if [ ! -f "$monolis_dir/src/matrix/sparse_util.f90" ]; then
        print_error "monolis source was not initialized correctly"
        exit 1
    fi

    if grep -q "subroutine monolis_add_scalar_to_sparse_matrix_atomic" "$monolis_dir/src/matrix/sparse_util.f90" \
        && grep -q -- "-fopenmp" "$monolis_dir/Makefile"; then
        print_success "monolis already contains the S-IGA atomic sparse-matrix addition and OpenMP flags"
        return
    fi

    if [ ! -f "$MONOLIS_PATCH" ]; then
        print_error "Missing monolis compatibility patch: $MONOLIS_PATCH"
        exit 1
    fi

    print_warning "monolis submodule lacks the required S-IGA compatibility changes. Applying project patch..."
    if git -C "$monolis_dir" apply --check "$MONOLIS_PATCH"; then
        git -C "$monolis_dir" apply "$MONOLIS_PATCH"
        print_success "Applied monolis S-IGA compatibility patch"
    elif git -C "$monolis_dir" apply --reverse --check "$MONOLIS_PATCH"; then
        print_success "monolis compatibility patch is already applied"
    else
        print_error "Failed to apply monolis compatibility patch"
        print_info "Long-term fix: push the monolis changes and update sfem_linear's submodule pointer."
        exit 1
    fi

    if grep -q "subroutine monolis_add_scalar_to_sparse_matrix_atomic" "$monolis_dir/src/matrix/sparse_util.f90" \
        && grep -q -- "-fopenmp" "$monolis_dir/Makefile"; then
        print_success "Verified monolis source compatibility changes"
    else
        print_error "monolis source compatibility changes are still missing after patch"
        exit 1
    fi
}

build_monolis() {
    local monolis_dir="$SOLVER_DIR/submodule/monolis"

    print_info "Building monolis dependencies..."
    pushd "$monolis_dir" > /dev/null
    ./install_lib.sh METIS

    print_info "Building monolis with MPI, METIS, and OpenMP support..."
    make clean > /dev/null 2>&1 || true
    make -B FLAGS=MPI,METIS

    if monolis_library_has_atomic_symbol "lib/libmonolis.a"; then
        print_success "Verified monolis atomic sparse-matrix symbol"
    else
        print_error "libmonolis.a does not contain the monolis atomic sparse-matrix module symbol"
        print_monolis_atomic_diagnostics
        popd > /dev/null
        exit 1
    fi
    popd > /dev/null
}

build_solver() {
    print_info "Building sfem_linear solver..."
    pushd "$SOLVER_DIR" > /dev/null
    make clean > /dev/null 2>&1 || true
    if make; then
        print_success "Solver built successfully!"
    else
        print_error "Failed to build solver!"
        popd > /dev/null
        exit 1
    fi
    popd > /dev/null
}

monolis_library_has_atomic_symbol() {
    local monolis_lib="${1:-$SOLVER_DIR/submodule/monolis/lib/libmonolis.a}"

    [ -f "$monolis_lib" ] && \
        nm -a "$monolis_lib" 2>/dev/null | grep -Eiq "monolis_add_scalar_to_sparse_matrix_atomic"
}

print_monolis_atomic_diagnostics() {
    print_info "Diagnostics for monolis atomic symbol:"
    print_info "mpif90 wrapper:"
    mpif90 --version | sed -n '1,3p' || true
    mpif90 -show 2>/dev/null || true

    print_info "Source check:"
    grep -n "monolis_add_scalar_to_sparse_matrix_atomic" src/matrix/sparse_util.f90 || true

    print_info "Object symbols containing atomic/add_scalar:"
    nm -a obj/matrix/sparse_util.o 2>/dev/null | grep -Ei "atomic|add_scalar" || true

    print_info "Archive symbols containing atomic/add_scalar:"
    nm -a lib/libmonolis.a 2>/dev/null | grep -Ei "atomic|add_scalar" || true
}

existing_monolis_library_has_atomic_symbol() {
    local monolis_lib="$SOLVER_DIR/submodule/monolis/lib/libmonolis.a"

    monolis_library_has_atomic_symbol "$monolis_lib"
}

prepare_sfem_linear_repo

REBUILD_SOLVER=0
if [ ! -f "$SOLVER_BIN" ]; then
    print_warning "Solver binary not found. Building solver and monolis..."
    REBUILD_SOLVER=1
elif [ "${FORCE_REBUILD_SOLVER:-0}" = "1" ]; then
    print_warning "FORCE_REBUILD_SOLVER=1; rebuilding solver and monolis"
    REBUILD_SOLVER=1
elif [ -f "$SOLVER_DIR/submodule/monolis/lib/libmonolis.a" ] && ! existing_monolis_library_has_atomic_symbol; then
    print_warning "Existing libmonolis.a lacks the atomic sparse-matrix symbol. Rebuilding solver and monolis..."
    REBUILD_SOLVER=1
fi

if [ "$REBUILD_SOLVER" -eq 1 ]; then
    if [ -f "$SOLVER_BIN" ] && [ "${FORCE_REBUILD_SOLVER:-0}" != "1" ]; then
        print_info "Rebuild is needed to keep future sfem_linear builds linkable"
    fi
    ensure_monolis_source
    build_monolis
    build_solver
else
    print_success "Fortran solver binary exists: $SOLVER_BIN"
fi

# Verify the binary is executable
if [ -f "$SOLVER_BIN" ]; then
    if [ -x "$SOLVER_BIN" ]; then
        print_success "Solver is ready to use"
    else
        print_warning "Making solver executable..."
        chmod +x "$SOLVER_BIN"
        print_success "Solver permissions updated"
    fi
else
    print_error "Solver binary still not found after build attempt!"
    print_info "You may need to manually build sfem_linear:"
    print_info "  FORCE_REBUILD_SOLVER=1 ./setup.sh"
fi

echo ""

# ============================================================================
# Installation Complete
# ============================================================================
print_header "Installation Complete! 🎉"

echo ""
print_success "Setup completed successfully!"
echo ""
print_info "To use the solver:"
echo ""
echo "  1. Activate the virtual environment:"
echo "     ${GREEN}pipenv shell${NC}"
echo ""
echo "  2. Run the solver from circular_crack directory:"
echo "     ${GREEN}cd circular_crack${NC}"
echo "     ${GREEN}python3 main.py --help${NC}"
echo ""
echo "  Alternative (without entering shell):"
echo "     ${GREEN}cd circular_crack${NC}"
echo "     ${GREEN}pipenv run python3 main.py --help${NC}"
echo ""
print_info "Virtual environment location: .venv/"
print_info "Python executable: $PYTHON_EXEC"
echo ""

# Show installed packages
print_info "Installed Python packages:"
pipenv run pip list | grep -E "(numpy|scipy|logzero)" || true
echo ""

print_success "You're all set! Happy computing! 🚀"
echo ""
