#!/bin/bash
# Generate Sneddon interpolation data for static mode
# This needs to be run ONCE before using --static_only mode

echo "========================================"
echo "Generating Sneddon Interpolation Data"
echo "========================================"
echo ""
echo "This script will precompute Sneddon analytical solution"
echo "using Bessel integrals (like Mathematica SneddonApp.mx)"
echo ""
echo "Options:"
echo "  1) Test mode (120x40 grid, ~5 minutes)"
echo "  2) Full mode (600x200 grid, ~30-45 minutes)"
echo ""
read -p "Select mode [1/2]: " mode

cd "$(dirname "$0")"

if [ "$mode" == "1" ]; then
    echo ""
    echo "Running in TEST mode..."
    python3 sneddon_precompute.py --test
    echo ""
    echo "✓ Test data generated: sneddon_interpolation_test.npz"
    echo ""
    echo "Note: This is low resolution. For production, run full mode."
    
elif [ "$mode" == "2" ]; then
    echo ""
    echo "Running in FULL mode..."
    echo "This will take ~30-45 minutes. Please be patient..."
    echo ""
    python3 sneddon_precompute.py
    echo ""
    echo "✓ Full data generated: sneddon_interpolation.npz"
    echo ""
    echo "You can now use --static_only mode!"
    
else
    echo "Invalid selection. Exiting."
    exit 1
fi

echo ""
echo "========================================"
echo "Done!"
echo "========================================"
