import numpy as np
from const import simulation_params as sp, material_property as mp
from utils.logger import logger
from utils.step2str import step2str


def generate(step, REstart):
    """
    Generate input.dat file for SFEM solver
    Format: simple value + comment format (matching Mathematica output)
    
    Based on Mathematica's input file structure from reference
    """
    logger.info(f"Generating input.dat for step {step}")
    
    from const import const_local_mesh as clm, const_global_mesh as cgm
    
    # Solution type: 1=static, 2=dynamic
    # In static mode, always use solution_type=1
    if sp.static_mode:
        solution_type = 1
        logger.info(f"Static mode: solution_type=1 (static analysis)")
    else:
        solution_type = 1 if step == 0 else 2
    
    # Nonlinear geometry: 0=off, 1=on
    is_nlgeom = 0
    
    # Max Newton-Raphson steps
    max_nr_step = 5
    
    # Time step size: dt = hL / V
    dt = clm.hL / sp.V
    
    # Max time steps
    max_time_step = 1
    
    # Number of increments
    nincrement = sp.inc
    
    # Material properties
    young = mp.EE
    poisson = mp.Nu
    density = mp.Rho
    
    # HHT-alpha parameter. The solver computes beta and gamma from alpha.
    alpha = sp.Alpha
    
    # Rayleigh damping
    Rm = sp.Alpha_l  # Mass damping coefficient (alpha)
    
    # Calculate Rk (Beta_l) using formula: β_l = 2.57 * h * sqrt(ρ/E)
    # islocal: 0=use hG (global mesh), 1=use hL (local mesh)
    # For S-IGA with local refinement, we use hL
    islocal = 1
    h = sp.hG if islocal == 0 else clm.hL
    Rk = 2.57 * h * np.sqrt(density / young)  # Stiffness damping coefficient (beta)
    
    # Local mesh parameters
    ngp = sp.ngp  # Integration points
    nrefine = sp.nrefLlist  # h-refine level
    
    # Restart flag
    # In static mode, always use is_restart=0 (no restart)
    if sp.static_mode:
        is_restart = 0
        logger.info(f"Static mode: is_restart=0 (fresh start)")
    else:
        is_restart = REstart
    
    # Number of threads
    num_threads = sp.OPENMP
    
    # IGA parameters
    order = cgm.p  # NURBS order
    nPtsU = sp.nPtsX  # Control points in U direction
    nPtsV = sp.nPtsY  # Control points in V direction
    nPtsW = sp.nPtsZ  # Control points in W direction
    
    # Penalty method flag
    penalty = 0
    is_sphere = 0
    
    # Write input.dat file (matching Mathematica format exactly)
    with open('input.dat', 'w') as f:
        # Line 1: solution type (1=static, 2=dynamic)
        # Note: Reference uses 1, but our analysis is dynamic. Use 2 for consistency with code
        f.write(f"{solution_type}\t  !>solutiontype(1:static 2:dynamic)\n")
        
        # Line 2: nonlinear geometry flag
        f.write(f"{is_nlgeom}\t  !>isNLgeom(0:off 1:on)\n")
        
        # Line 3: max NR steps
        f.write(f"{max_nr_step}\t  !> max NR step\n")
        
        # Line 4: time step (use scientific notation like reference)
        f.write(f"{dt:.1e}\t  !> dt\n")
        
        # Line 5: max time step
        f.write(f"{max_time_step}\t!>max time step\n")
        
        # Line 6: number of increments
        f.write(f"{nincrement}\t  !> nincrement\n")
        
        # Line 7: Young's modulus
        f.write(f"{young:.2e}\t  !> young 率\n")
        
        # Line 8: Poisson ratio
        f.write(f"{poisson}\t  !> Poisson 比\n")
        
        # Line 9: density (with trailing dot like reference)
        f.write(f"{density:.0f}.\t  !> density\n")
        
        # Line 10: HHT alpha
        f.write(f"{alpha}\t  !> alpha for HHT-alpha\n")
        
        # Line 11: Rm (with trailing dot like reference)
        f.write(f"{Rm:.0f}.\t  !> Rm for Rayleigh damping for mass\n")
        
        # Line 12: Rk (full precision scientific notation)
        f.write(f"{Rk}\t  !> Rk for Rayleigh damping for stiffness\n")
        
        # Line 13: integration points
        f.write(f"{ngp}\t  !> local mesh integral point (not used)\n")
        
        # Line 14: h-refine
        f.write(f"{nrefine}\t  !> local mesh h-refine\n")
        
        # Line 15: restart flag
        f.write(f"{is_restart}\t  !> is_Restart (0:off, 1:on)\n")
        
        # Line 16: number of threads
        f.write(f"{num_threads}\t  !> number of thread\n")
        
        # Line 17: IGA order
        f.write(f"{order}\t  !> order for IGA\n")
        
        # Line 18: control points U
        f.write(f"{nPtsU}\t  !> number of controlPts in U-direction\n")
        
        # Line 19: control points V
        f.write(f"{nPtsV}\t  !> number of controlPts in V-direction\n")
        
        # Line 20: control points W
        f.write(f"{nPtsW}\t  !> number of controlPts in W-direction\n")
        
        # Line 21: penalty flag (always write, value=0)
        f.write(f"{penalty}\t  !> penalty or not\n")

        # Line 22: is_sphere flag (always write, value=0)
        f.write(f"{is_sphere}\t  !> is_sphere or not\n")
    
    logger.info("input.dat generated successfully")


def generate_virtual_mesh():
    """
    Virtual mesh files are generated by GlobalMesh.generate()
    This function is kept for compatibility but does nothing
    """
    # node.v.dat and elem.v.dat are already generated by GlobalMesh._build_visual_mesh()
    # Do not overwrite them here
    logger.info("Virtual mesh files already generated by GlobalMesh")
