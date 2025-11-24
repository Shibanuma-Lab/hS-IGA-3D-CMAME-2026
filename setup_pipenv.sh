#!/bin/bash
# Setup script using pipenv

echo "=========================================="
echo "S-IGA Setup with Pipenv"
echo "=========================================="
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

# Check if pipenv is installed
echo "Checking for pipenv..."
if ! command -v pipenv &> /dev/null; then
    echo "❌ pipenv not found. Installing..."
    
    # Try to install pipenv
    if command -v pip3 &> /dev/null; then
        pip3 install --user pipenv
    elif command -v pip &> /dev/null; then
        pip install --user pipenv
    else
        echo "❌ Error: pip not found. Install pip first:"
        echo "  sudo apt install python3-pip"
        echo "  pip3 install --user pipenv"
        exit 1
    fi
    
    # Add to PATH if needed
    if ! command -v pipenv &> /dev/null; then
        echo ""
        echo "⚠️  pipenv installed but not in PATH"
        echo "Add this to your ~/.bashrc:"
        echo '  export PATH="$HOME/.local/bin:$PATH"'
        echo ""
        echo "Then run: source ~/.bashrc"
        exit 1
    fi
fi
echo "✅ pipenv found"
echo ""

# Install dependencies
echo "Installing dependencies with pipenv..."
pipenv install
if [ $? -ne 0 ]; then
    echo "❌ Error: Failed to install dependencies"
    exit 1
fi
echo "✅ Dependencies installed"
echo ""

# Verify installation
echo "Verifying installation..."
pipenv run python -c "import numpy, scipy, logzero; print('✅ All packages imported successfully')"
if [ $? -ne 0 ]; then
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
echo "To use the simulation:"
echo "  1. Enter pipenv shell:"
echo "     pipenv shell"
echo ""
echo "  2. Run the simulation:"
echo "     cd circular_crack"
echo "     python main.py --help"
echo ""
echo "  3. Exit pipenv shell:"
echo "     exit"
echo ""
echo "Or run directly without shell:"
echo "  pipenv run python circular_crack/main.py --help"
echo ""
