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
# ============================================================================

set -e  # Exit on error

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
    
    print_info "Installing build dependencies..."
    sudo apt update
    sudo apt upgrade -y
    sudo apt install -y build-essential zlib1g-dev libncurses5-dev libgdbm-dev \
                        libnss3-dev libssl-dev libsqlite3-dev libreadline-dev \
                        libffi-dev libbz2-dev wget
    
    print_info "Downloading Python $PYTHON_VERSION..."
    cd /tmp
    wget https://www.python.org/ftp/python/$PYTHON_VERSION/Python-$PYTHON_VERSION.tgz
    
    print_info "Extracting and building Python..."
    tar -xf Python-$PYTHON_VERSION.tgz
    cd Python-$PYTHON_VERSION
    ./configure --enable-optimizations
    make -j$(nproc)
    
    print_info "Installing Python 3.10..."
    sudo make altinstall
    
    # Clean up
    cd /tmp
    rm -rf Python-$PYTHON_VERSION Python-$PYTHON_VERSION.tgz
    
    # Verify installation
    if python3.10 --version &> /dev/null; then
        print_success "Python 3.10 installed successfully!"
        PYTHON_EXEC=$(command -v python3.10)
    else
        print_error "Python 3.10 installation failed!"
        exit 1
    fi
    
    # Return to project directory
    cd "$OLDPWD"
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
    
    print_info "Installing Pipenv..."
    python3.10 -m pip install --user pipenv
    
    # Add to PATH in .bashrc if not already present
    BASHRC="$HOME/.bashrc"
    PATH_EXPORT='export PATH="$HOME/.local/bin:$PATH"'
    
    if ! grep -q "$PATH_EXPORT" "$BASHRC"; then
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
pipenv install

if [ $? -eq 0 ]; then
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

SOLVER_BIN="sfem_linear/bin/sfem_linear"

# Check if sfem_linear directory exists
if [ ! -d "sfem_linear" ]; then
    print_warning "sfem_linear directory not found. Cloning from GitLab..."
    
    # Check if git is available
    if ! command -v git &> /dev/null; then
        print_error "git is not installed!"
        print_info "Please install git: sudo apt install git"
        exit 1
    fi
    
    print_info "Cloning sfem_linear repository..."
    if git clone git@gitlab.com:morita/sfem_linear.git; then
        print_success "sfem_linear repository cloned successfully"
    else
        print_error "Failed to clone sfem_linear repository!"
        print_info "Please ensure you have SSH access to gitlab.com:morita/sfem_linear.git"
        print_info "Or manually clone the repository and re-run this script"
        exit 1
    fi
    
    # Enter the directory and checkout the correct branch
    cd sfem_linear
    print_info "Switching to branch: tianyu_IGA"
    if git checkout tianyu_IGA; then
        print_success "Switched to branch tianyu_IGA"
    else
        print_error "Failed to checkout branch tianyu_IGA!"
        cd ..
        exit 1
    fi
    
    # Run install script if it exists
    if [ -f "install_lib.sh" ]; then
        print_info "Running install_lib.sh..."
        if bash install_lib.sh; then
            print_success "install_lib.sh completed"
        else
            print_warning "install_lib.sh encountered issues (this may be normal)"
        fi
    else
        print_warning "install_lib.sh not found, skipping..."
    fi
    
    # Build the solver (sequential build to avoid Fortran module dependency issues)
    print_info "Building solver with make..."
    if make; then
      print_success "Solver built successfully!"
    else
      print_error "Failed to build solver!"
      print_info "Please check the build errors above"
      cd ..
      exit 1
    fi    # Return to project root
    cd ..
    
elif [ ! -f "$SOLVER_BIN" ]; then
    # Directory exists but binary doesn't - try to build
    print_warning "Solver binary not found. Attempting to build..."
    
    cd sfem_linear
    
    # Check if we're on the correct branch
    CURRENT_BRANCH=$(git branch --show-current 2>/dev/null || echo "unknown")
    if [ "$CURRENT_BRANCH" != "tianyu_IGA" ]; then
        print_info "Switching to branch: tianyu_IGA"
        git checkout tianyu_IGA || print_warning "Could not switch branch"
    fi
    
    # Run install script if it exists
    if [ -f "install_lib.sh" ]; then
        print_info "Running install_lib.sh..."
        bash install_lib.sh || print_warning "install_lib.sh encountered issues"
    fi
    
    # Build the solver (sequential build to avoid Fortran module dependency issues)
    print_info "Building solver with make..."
    if make; then
      print_success "Solver built successfully!"
    else
      print_error "Failed to build solver!"
      cd ..
      exit 1
    fi
    
    cd ..
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
    print_info "  cd sfem_linear"
    print_info "  git checkout tianyu_IGA"
    print_info "  bash install_lib.sh"
    print_info "  make"
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
