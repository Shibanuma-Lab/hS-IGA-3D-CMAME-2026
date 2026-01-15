#!/bin/bash
# Script to ignore local changes to const files
# These files are tracked in git but their local modifications will be ignored

echo "Ignoring local changes to const files..."

# List of const files to ignore
const_files=(
    "circular_crack/const/const_global_mesh.py"
    "circular_crack/const/const_local_mesh.py"
    "circular_crack/const/const_jintegral.py"
    "circular_crack/const/simulation_params.py"
    "circular_crack/const/material_property.py"
)

for file in "${const_files[@]}"; do
    if [ -f "$file" ]; then
        git update-index --skip-worktree "$file"
        echo "✓ Ignoring changes to: $file"
    else
        echo "⚠ File not found: $file"
    fi
done

echo ""
echo "Done! Local changes to these files will now be ignored."
echo "To see which files are being ignored, run: git ls-files -v | grep '^S'"
