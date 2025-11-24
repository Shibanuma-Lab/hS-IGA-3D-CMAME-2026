#!/bin/bash
# Quick setup script for S-IGA Circular Crack Simulation

echo "=========================================="
echo "S-IGA Setup Script"
echo "=========================================="
echo ""

# Check Python version
echo "Checking Python version..."
python3 --version
if [ $? -ne 0 ]; then
    echo "❌ Error: Python 3 is not installed"
    exit 1
fi
echo "✅ Python 3 found"
echo ""

# Install Python dependencies
echo "Installing Python dependencies..."
pip3 install -r requirements.txt
if [ $? -ne 0 ]; then
    echo "❌ Error: Failed to install Python dependencies"
    echo "Try running: pip3 install --user -r requirements.txt"
    exit 1
fi
echo "✅ Python dependencies installed"
echo ""

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
