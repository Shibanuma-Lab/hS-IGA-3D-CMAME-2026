#!/bin/bash
# Quick setup script for S-IGA Circular Crack Simulation

echo "=========================================="
echo "S-IGA Setup Script"
echo "=========================================="
echo ""

# Detect WSL environment
if grep -qi microsoft /proc/version; then
    echo "🐧 Detected WSL (Windows Subsystem for Linux)"
    WSL_ENV=true
else
    WSL_ENV=false
fi
echo ""

# Check Python version
echo "Checking Python version..."
python3 --version
if [ $? -ne 0 ]; then
    echo "❌ Error: Python 3 is not installed"
    echo ""
    echo "Install Python 3 with:"
    echo "  sudo apt update"
    echo "  sudo apt install python3"
    exit 1
fi
echo "✅ Python 3 found"
echo ""

# Install Python dependencies
echo "Installing Python dependencies..."

# Try different pip commands and verify they actually work
PIP_CMD=""
if command -v pip3 &> /dev/null; then
    pip3 --version &> /dev/null && PIP_CMD="pip3"
fi

if [ -z "$PIP_CMD" ] && command -v pip &> /dev/null; then
    pip --version &> /dev/null && PIP_CMD="pip"
fi

if [ -z "$PIP_CMD" ]; then
    python3 -m pip --version &> /dev/null && PIP_CMD="python3 -m pip"
fi

if [ -z "$PIP_CMD" ]; then
    echo "❌ Error: pip is not installed"
    echo ""
    if [ "$WSL_ENV" = true ]; then
        echo "🔧 Quick fix for WSL:"
        echo ""
        echo "Option 1 - Install via apt (Recommended, fastest):"
        echo "  sudo apt update"
        echo "  sudo apt install python3-numpy python3-scipy"
        echo ""
        echo "Option 2 - Install pip first, then use pip:"
        echo "  sudo apt update"
        echo "  sudo apt install python3-pip"
        echo "  pip3 install --user numpy scipy"
    else
        echo "Please install pip first:"
        echo "  Ubuntu/Debian: sudo apt install python3-pip"
        echo "  macOS: brew install python3"
        echo ""
        echo "Or install packages manually:"
        echo "  python3 -m ensurepip --default-pip"
        echo "  pip3 install numpy scipy"
    fi
    exit 1
fi

echo "Using pip command: $PIP_CMD"

# Check if packages are already installed
echo "Checking existing packages..."
python3 -c "import numpy, scipy" 2>/dev/null
if [ $? -eq 0 ]; then
    echo "✅ NumPy and SciPy already installed"
    echo ""
else
    echo "Installing packages..."
    $PIP_CMD install -r requirements.txt
    if [ $? -ne 0 ]; then
        echo "⚠️  System-wide installation failed. Trying user installation..."
        $PIP_CMD install --user -r requirements.txt
        if [ $? -ne 0 ]; then
            echo "❌ Error: Failed to install Python dependencies"
            echo ""
            if [ "$WSL_ENV" = true ]; then
                echo "WSL-specific troubleshooting:"
                echo "  1. Update package list:"
                echo "     sudo apt update"
                echo "  2. Install pip3:"
                echo "     sudo apt install python3-pip"
                echo "  3. Install packages:"
                echo "     pip3 install --user numpy scipy"
                echo ""
                echo "  Alternative: Use apt packages (faster in WSL):"
                echo "     sudo apt install python3-numpy python3-scipy"
            else
                echo "Please try manually:"
                echo "  sudo apt-get install python3-pip  # Install pip first"
                echo "  pip3 install --user numpy scipy   # Install packages"
            fi
            exit 1
        fi
    fi
    echo "✅ Python dependencies installed"
    echo ""
    
    # Verify installation
    echo "Verifying installation..."
    python3 -c "import numpy; import scipy; print(f'  NumPy version: {numpy.__version__}'); print(f'  SciPy version: {scipy.__version__}')"
    if [ $? -eq 0 ]; then
        echo "✅ All packages working correctly"
    else
        echo "⚠️  Warning: Packages installed but import failed"
        echo "   You may need to restart your terminal or run:"
        echo "   export PATH=\"\$HOME/.local/bin:\$PATH\""
    fi
    echo ""
fi

# Check for Fortran compiler
echo "Checking for Fortran compiler..."
which gfortran > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "⚠️  Warning: gfortran not found. Please install a Fortran compiler."
    echo "   Ubuntu/Debian: sudo apt-get install gfortran"
    echo "   macOS: brew install gcc"
else
    echo "✅ Fortran compiler found"
fi
echo ""

# Build solver
if [ -d "sfem_linear" ]; then
    echo "Building Fortran solver..."
    cd sfem_linear
    make
    if [ $? -ne 0 ]; then
        echo "❌ Error: Failed to build solver"
        cd ..
        exit 1
    fi
    cd ..
    echo "✅ Solver built successfully"
    echo ""
    
    # Verify solver
    if [ -f "sfem_linear/bin/sfem_linear" ]; then
        echo "✅ Solver executable found at: sfem_linear/bin/sfem_linear"
    else
        echo "⚠️  Warning: Solver executable not found"
    fi
else
    echo "⚠️  Warning: sfem_linear directory not found. Skipping solver build."
fi

echo ""
echo "=========================================="
echo "Setup Complete! 🎉"
echo "=========================================="
echo ""
echo "To run the simulation:"
echo "  cd circular_crack"
echo "  python3 main.py --help"
echo ""
