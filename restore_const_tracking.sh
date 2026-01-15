#!/bin/bash
# Script to restore tracking of const files
# Use this when you want to commit changes to const files

echo "Restoring tracking of const files..."

# List of const files to restore
const_files=(
    "circular_crack/const/const_global_mesh.py"
    "circular_crack/const/const_local_mesh.py"
    "circular_crack/const/const_jintegral.py"
    "circular_crack/const/simulation_params.py"
    "circular_crack/const/material_property.py"
)

for file in "${const_files[@]}"; do
    if [ -f "$file" ]; then
        git update-index --no-skip-worktree "$file"
        echo "✓ Restored tracking: $file"
    else
        echo "⚠ File not found: $file"
    fi
done

echo ""
echo "Done! Changes to these files will now be tracked again."
echo "You can now commit changes if needed."
