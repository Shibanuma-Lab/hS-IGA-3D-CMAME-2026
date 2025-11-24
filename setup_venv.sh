#!/bin/bash
# Setup script using virtual environment (venv)

echo "=========================================="
echo "S-IGA Setup with Virtual Environment"
echo "=========================================="
echo ""

# Detect WSL
if grep -qi microsoft /proc/version; then
    echo "🐧 Detected WSL (Windows Subsystem for Linux)"
    WSL_ENV=true
else
    WSL_ENV=false
fi
echo ""

# Check Python
echo "Checking Python version..."
python3 --version
if [ $? -ne 0 ]; then
    echo "❌ Error: Python 3 is not installed"
    exit 1
fi
echo "✅ Python 3 found"
echo ""

# Check if venv module is available
echo "Checking for venv module..."
python3 -m venv --help > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "❌ Error: venv module not found"
    echo ""
    echo "Install it with:"
    echo "  sudo apt install python3-venv"
    exit 1
fi
echo "✅ venv module found"
echo ""

# Create virtual environment
if [ -d "venv" ]; then
    echo "⚠️  Virtual environment already exists at ./venv"
    read -p "Do you want to recreate it? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf venv
        echo "Removed old virtual environment"
    else
        echo "Using existing virtual environment"
    fi
fi

if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
    if [ $? -ne 0 ]; then
        echo "❌ Error: Failed to create virtual environment"
        exit 1
    fi
    echo "✅ Virtual environment created at ./venv"
fi
echo ""

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate
echo "✅ Virtual environment activated"
echo ""

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip
echo ""

# Install dependencies
echo "Installing Python dependencies..."
pip install -r requirements.txt
if [ $? -ne 0 ]; then
    echo "❌ Error: Failed to install dependencies"
    exit 1
fi
echo "✅ Dependencies installed"
echo ""

# Verify installation
echo "Verifying installation..."
python -c "import numpy, scipy, logzero; print(f'  NumPy: {numpy.__version__}'); print(f'  SciPy: {scipy.__version__}'); print(f'  logzero: {logzero.__version__}')"
if [ $? -eq 0 ]; then
    echo "✅ All packages working correctly"
else
    echo "❌ Error: Package verification failed"
    exit 1
fi
echo ""

# Check Fortran compiler
echo "Checking for Fortran compiler..."
which gfortran > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "⚠️  Warning: gfortran not found"
    echo "   Install with: sudo apt install gfortran"
else
    echo "✅ Fortran compiler found"
    
    # Build solver
    if [ -d "sfem_linear" ]; then
        echo ""
        echo "Building Fortran solver..."
        cd sfem_linear
        make
        if [ $? -eq 0 ]; then
            echo "✅ Solver built successfully"
            cd ..
        else
            echo "❌ Error: Failed to build solver"
            cd ..
            exit 1
        fi
    fi
fi
echo ""

echo "=========================================="
echo "Setup Complete! 🎉"
echo "=========================================="
echo ""
echo "Virtual environment is located at: ./venv"
echo ""
echo "To use the simulation:"
echo "  1. Activate virtual environment:"
echo "     source venv/bin/activate"
echo ""
echo "  2. Run the simulation:"
echo "     cd circular_crack"
echo "     python main.py --help"
echo ""
echo "  3. When done, deactivate:"
echo "     deactivate"
echo ""

# Deactivate for now
deactivate
