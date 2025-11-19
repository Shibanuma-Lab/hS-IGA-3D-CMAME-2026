"""
Workflow Diagram for Circular Crack Simulation
===============================================

                         START
                           |
                           v
                   +--------------+
                   | Configuration|
                   | (const/*.py) |
                   +--------------+
                           |
                           v
        +------------------+------------------+
        |                                     |
        v                                     v
   [Step Loop: 0 to stepall-1]          [First Run?]
        |                                     |
        |                                     v
        v                                  [Setup]
   +--------------------+              - Create directories
   | Generate Meshes    |              - Initialize logger
   +--------------------+              - Check solver exists
        |
        |---> GlobalMesh.make_global_mesh()
        |     - Graded refinement (m0 levels)
        |     - Crack length: a = c*step/stepall
        |     - Output: node.g.dat, elem.g.dat, weights.g.dat
        |
        |---> LocalMesh.make_local_mesh()
        |     - Polar coordinates around crack
        |     - Fine elements near tip
        |     - Output: node.l.dat, elem.l.dat, weights.l.dat
        |
        v
   +--------------------+
   | Define Boundary    |
   +--------------------+
        |
        |---> Boundary.define_boundary()
        |     - Symmetry planes (ux=0, uy=0, uz=0)
        |     - Remote stress as nodal forces
        |     - Output: bc.g.dat, bc.l.dat, load.dat
        |
        v
   +--------------------+
   | Generate Input     |
   +--------------------+
        |
        |---> input_generator.generate()
        |     - Material properties
        |     - Time integration (HHT)
        |     - Mesh file references
        |     - Output: input.dat
        |
        |---> input_generator.generate_virtual_mesh()
              - Output: node.v.dat, elem.v.dat
        |
        v
   +--------------------+
   | Check REstart?     |
   +--------------------+
        |
        +---> [step > 0] --> initial()
        |                    - Read delta_u.dat, velocity.dat
        |                    - from previous step
        |                    - Initialize current step
        |
        +---> [step == 0] -> Zero initial conditions
        |
        v
   +--------------------+
   | Run Solver         |
   +--------------------+
        |
        |---> linux_command.run()
        |     - Execute: sfem_linear input.dat
        |     - Set OMP_NUM_THREADS
        |     - Capture output
        |
        v
   +--------------------+
   | Check Results      |
   +--------------------+
        |
        |---> Check if output files exist:
        |     - delta_u.dat (displacement)
        |     - velocity.dat
        |     - acceleration.dat
        |     - stress.dat
        |     - strain.dat
        |
        v
   +--------------------+
   | Extract Results    |
   +--------------------+
        |
        |---> Copy to output folder
        |---> Prepare for next step
        |
        v
   +--------------------+
   | Next Step?         |
   +--------------------+
        |
        +---> [step < stepall-1] --+
        |                          |
        v                          |
      [END]  <--------------------+


Data Flow Diagram
=================

Configuration Files:
  const_global_mesh.py --+
  const_local_mesh.py ---+---> Parameters
  material_property.py --+
  simulation_params.py --+
                         |
                         v
                   +----------+
                   |  main.py |
                   +----------+
                         |
        +----------------+----------------+
        |                |                |
        v                v                v
  GlobalMesh       LocalMesh        Boundary
        |                |                |
        v                v                v
   node.g.dat       node.l.dat        bc.g.dat
   elem.g.dat       elem.l.dat        bc.l.dat
   weights.g.dat    weights.l.dat     load.dat
        |                |                |
        +----------------+----------------+
                         |
                         v
                input_generator
                         |
                         v
                    input.dat
                         |
                         v
               +------------------+
               | sfem_linear      |
               | (Fortran Solver) |
               +------------------+
                         |
        +----------------+----------------+
        |                |                |
        v                v                v
   delta_u.dat    velocity.dat      stress.dat
   acceleration.dat  strain.dat   reaction.dat


File Format Examples
====================

node.g.dat (Global Mesh Nodes):
  NodeID    X-coord       Y-coord       Z-coord
  -------   -----------   -----------   -----------
  1         0.000000e+00  -5.000000e-01  0.000000e+00
  2         5.000000e-04  -5.000000e-01  0.000000e+00
  ...

elem.g.dat (Global Mesh Elements):
  ElemID  Node1  Node2  Node3  ...  Node64
  ------  -----  -----  -----  ---  ------
  1       1      2      3      ...  64
  2       2      3      4      ...  65
  ...

bc.g.dat (Boundary Conditions):
  NodeID  DOF   Value
  ------  ---   -----
  1       1     0.0      # ux = 0
  15      2     0.0      # uy = 0
  ...

load.dat (Applied Loads):
  NodeID  DOF   Force
  ------  ---   -------
  100     1     1.5e6   # Fx
  101     1     1.5e6
  ...

input.dat (Solver Input):
  # Problem type
  DYNAMIC
  
  # Mesh files
  GLOBAL_MESH
    NODE_FILE      node.g.dat
    ELEMENT_FILE   elem.g.dat
  ...
  
  # Material properties
  MATERIAL
    YOUNGS_MODULUS   2.06e11
    POISSON_RATIO    0.3
  ...


Parameter Sensitivity
=====================

Critical Parameters:
  rGL (Global/Local ratio):
    - Too small: Excessive DOF
    - Too large: Poor resolution
    - Recommended: 4-10

  hL (Local element size):
    - Should resolve crack tip field
    - Typical: c/100 to c/50
    - Trade-off: accuracy vs. cost

  m0 (Refinement levels):
    - More levels: Smoother transition
    - Fewer levels: Fewer DOF
    - Recommended: 3-5

  OPENMP (Threads):
    - Set to physical cores
    - Avoid hyperthreading
    - Monitor memory bandwidth


Typical Simulation Timeline
============================

Step 0 (Initial):
  ├─ Generate meshes           (< 1 sec)
  ├─ Define BCs                (< 1 sec)
  ├─ Create input files        (< 1 sec)
  └─ Run solver                (minutes to hours)
       └─ Assembly & solve
       └─ Write results

Step 1-N (Propagation):
  ├─ Generate new meshes       (< 1 sec)
  │    └─ Crack length updated
  ├─ Define BCs                (< 1 sec)
  ├─ Create input files        (< 1 sec)
  ├─ Load previous results     (< 1 sec)
  └─ Run solver                (minutes to hours)
       └─ Read restart data
       └─ Time integration
       └─ Write results


Memory Requirements
===================

Global Mesh:
  Nodes:    nx × ny × nz
  Elements: (nx-p) × (ny-p) × (nz-p)
  Memory:   ~8 bytes × nnodes × 3 coords
  
Local Mesh:
  Nodes:    nr × nθ × nz
  Elements: (nr-p) × nθ × (nz-p)
  Memory:   Similar to global

Typical Case (10,000 nodes):
  Node storage:     ~240 KB
  Element storage:  ~1 MB
  BC storage:       ~100 KB
  Total input:      ~2 MB per step

Solver Memory:
  Stiffness matrix: ~sparse storage
  Working arrays:   ~DOF × 10
  Typical:          1-10 GB for 100K DOF


Error Checking
==============

Common Issues:
  
  1. "SFEM executable not found"
     → Compile solver: cd ../sfem_linear && make
  
  2. "Missing output files"
     → Check solver log
     → Verify input.dat format
  
  3. "Mesh generation failed"
     → Check parameters (positive values)
     → Verify crack length < domain size
  
  4. "Singular matrix"
     → Check boundary conditions
     → Ensure sufficient constraints
  
  5. "Excessive memory"
     → Reduce mesh density
     → Increase refinement ratio (rGL)

"""

print(__doc__)
