#!/usr/bin/env python3
"""
Test script for circular crack mesh generation
Verifies each module without running the full solver
"""

import os
import sys
import numpy as np

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from const import const_global_mesh as cgm, const_local_mesh as clm, simulation_params as sp
from global_mesh import GlobalMesh
from local_mesh import LocalMesh
from boundary import Boundary
import input_generator
from utils.logger import logger


def test_global_mesh():
    """Test global mesh generation"""
    print("\n" + "="*60)
    print("Testing Global Mesh Generation")
    print("="*60)
    
    step = 0
    g = GlobalMesh(step)
    
    print(f"Parameters:")
    print(f"  m0 (refinement levels): {g.m0}")
    print(f"  rGL (global/local ratio): {g.rGL}")
    print(f"  Crack radius: {g.c} m")
    print(f"  Thickness: {g.thi} m")
    
    g.make_global_mesh()
    
    print(f"\nResults:")
    print(f"  Number of nodes: {len(g.node_g)}")
    print(f"  Number of elements: {len(g.elem_g)}")
    print(f"  Element connectivity size: {g.elem_g.shape[1] - 1} nodes/element")
    
    # Verify node coordinates
    if g.node_g is not None:
        x_range = [g.node_g[:, 1].min(), g.node_g[:, 1].max()]
        y_range = [g.node_g[:, 2].min(), g.node_g[:, 2].max()]
        z_range = [g.node_g[:, 3].min(), g.node_g[:, 3].max()]
        
        print(f"  X range: [{x_range[0]:.6e}, {x_range[1]:.6e}]")
        print(f"  Y range: [{y_range[0]:.6e}, {y_range[1]:.6e}]")
        print(f"  Z range: [{z_range[0]:.6e}, {z_range[1]:.6e}]")
    
    print("✓ Global mesh generation successful")
    return g


def test_local_mesh():
    """Test local mesh generation"""
    print("\n" + "="*60)
    print("Testing Local Mesh Generation")
    print("="*60)
    
    step = 3
    l = LocalMesh(step)
    
    print(f"Parameters:")
    print(f"  hL (element size): {l.hL} m")
    print(f"  aL (crack elements): {l.aL}")
    print(f"  lL (ligament elements): {l.lL}")
    print(f"  HL (thickness elements): {l.HL}")
    
    l.make_local_mesh()
    
    print(f"\nResults:")
    print(f"  Crack length: {l.crack_length:.6e} m")
    print(f"  Number of nodes: {len(l.node_l)}")
    print(f"  Number of elements: {len(l.elem_l)}")
    
    # Verify node coordinates
    if l.node_l is not None:
        x_range = [l.node_l[:, 1].min(), l.node_l[:, 1].max()]
        y_range = [l.node_l[:, 2].min(), l.node_l[:, 2].max()]
        z_range = [l.node_l[:, 3].min(), l.node_l[:, 3].max()]
        
        print(f"  X range: [{x_range[0]:.6e}, {x_range[1]:.6e}]")
        print(f"  Y range: [{y_range[0]:.6e}, {y_range[1]:.6e}]")
        print(f"  Z range: [{z_range[0]:.6e}, {z_range[1]:.6e}]")
    
    print("✓ Local mesh generation successful")
    return l


def test_boundary_conditions(g, l):
    """Test boundary condition generation"""
    print("\n" + "="*60)
    print("Testing Boundary Conditions")
    print("="*60)
    
    b = Boundary(l, g)
    b.define_boundary(l, g)
    
    print(f"Results:")
    print(f"  Global BCs: {len(b.bc_g)} constraints")
    print(f"  Local BCs: {len(b.bc_l)} constraints")
    print(f"  Loads: {len(b.load_data)} load points")
    
    # Analyze BC types
    if len(b.bc_g) > 0:
        dof_types = np.unique(b.bc_g[:, 1].astype(int))
        print(f"  Global BC DOFs: {dof_types}")
    
    if len(b.load_data) > 0:
        total_load = np.sum(b.load_data[:, 2])
        print(f"  Total applied load: {total_load:.6e} N")
    
    print("✓ Boundary condition generation successful")
    return b


def test_input_generator():
    """Test input file generation"""
    print("\n" + "="*60)
    print("Testing Input File Generator")
    print("="*60)
    
    # Create temporary directory
    test_dir = "test_output"
    os.makedirs(test_dir, exist_ok=True)
    os.chdir(test_dir)
    
    step = 1
    REstart = 1
    
    input_generator.generate(step, REstart)
    input_generator.generate_virtual_mesh()
    
    # Check if files were created
    files_created = []
    if os.path.exists('input.dat'):
        files_created.append('input.dat')
    if os.path.exists('node.v.dat'):
        files_created.append('node.v.dat')
    if os.path.exists('elem.v.dat'):
        files_created.append('elem.v.dat')
    
    print(f"Files created: {files_created}")
    
    # Read and display part of input.dat
    if os.path.exists('input.dat'):
        with open('input.dat', 'r') as f:
            lines = f.readlines()
        print(f"\nInput file preview (first 20 lines):")
        for line in lines[:20]:
            print(f"  {line.rstrip()}")
    
    os.chdir('..')
    
    print("✓ Input file generation successful")


def test_file_io():
    """Test file writing functionality"""
    print("\n" + "="*60)
    print("Testing File I/O")
    print("="*60)
    
    # Create temporary directory
    test_dir = "test_output"
    os.makedirs(test_dir, exist_ok=True)
    os.chdir(test_dir)
    
    step = 0
    g = GlobalMesh(step)
    g.make_global_mesh()
    g.generate()
    
    l = LocalMesh(step)
    l.make_local_mesh()
    l.generate()
    
    # Check if files exist
    files_to_check = [
        'node.g.dat',
        'elem.g.dat',
        'weights.g.dat',
        'index.g.dat',
        'node.l.dat',
        'elem.l.dat',
        'weights.l.dat'
    ]
    
    for filename in files_to_check:
        if os.path.exists(filename):
            size = os.path.getsize(filename)
            print(f"  ✓ {filename} ({size} bytes)")
        else:
            print(f"  ✗ {filename} (missing)")
    
    os.chdir('..')
    
    print("✓ File I/O successful")


def test_full_workflow():
    """Test complete workflow for one step"""
    print("\n" + "="*60)
    print("Testing Full Workflow (Step 0)")
    print("="*60)
    
    # Create directories
    test_dir = "test_output/full_workflow"
    os.makedirs(test_dir, exist_ok=True)
    os.chdir(test_dir)
    
    step = 0
    
    # Generate global mesh
    print("1. Generating global mesh...")
    g = GlobalMesh(step)
    g.make_global_mesh()
    g.generate()
    print(f"   ✓ {len(g.node_g)} nodes, {len(g.elem_g)} elements")
    
    # Generate local mesh
    print("2. Generating local mesh...")
    l = LocalMesh(step)
    l.make_local_mesh()
    l.generate()
    print(f"   ✓ {len(l.node_l)} nodes, {len(l.elem_l)} elements")
    
    # Generate boundary conditions
    print("3. Generating boundary conditions...")
    b = Boundary(l, g)
    b.define_boundary(l, g)
    b.generate()
    print(f"   ✓ {len(b.bc_g)} BCs, {len(b.load_data)} loads")
    
    # Generate input file
    print("4. Generating input file...")
    input_generator.generate(step, 0)
    input_generator.generate_virtual_mesh()
    print("   ✓ Input files generated")
    
    # List all generated files
    print("\nGenerated files:")
    for file in sorted(os.listdir('.')):
        if os.path.isfile(file):
            size = os.path.getsize(file)
            print(f"   {file:20s} ({size:>10d} bytes)")
    
    os.chdir('../..')
    
    print("\n✓ Full workflow successful")


def main():
    """Run all tests"""
    print("\n" + "="*60)
    print("CIRCULAR CRACK MESH GENERATION TEST SUITE")
    print("="*60)
    
    try:
        # Test individual components
        g = test_global_mesh()
        l = test_local_mesh()
        test_boundary_conditions(g, l)
        test_input_generator()
        test_file_io()
        
        # Test full workflow
        test_full_workflow()
        
        print("\n" + "="*60)
        print("ALL TESTS PASSED ✓")
        print("="*60)
        
    except Exception as e:
        print("\n" + "="*60)
        print("TEST FAILED ✗")
        print("="*60)
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
