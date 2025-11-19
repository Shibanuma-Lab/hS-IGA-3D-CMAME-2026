import numpy as np
from const import simulation_params as sp, material_property as mp
from utils.logger import logger
from utils.step2str import step2str


def generate(step, REstart):
    """
    Generate input file for SFEM solver
    Writes all simulation parameters and settings
    """
    logger.info(f"Generating input.dat for step {step}")
    
    str_step = step2str(step)
    
    with open('input.dat', 'w') as f:
        # Write header
        f.write("# Input file for SFEM linear solver\n")
        f.write(f"# Step: {step}\n")
        f.write(f"# REstart: {REstart}\n")
        f.write("\n")
        
        # Problem type
        f.write("# Problem type\n")
        f.write("DYNAMIC\n")  # or STATIC for static analysis
        f.write("\n")
        
        # Mesh files
        f.write("# Mesh files\n")
        f.write("GLOBAL_MESH\n")
        f.write(f"  NODE_FILE      node.g.dat\n")
        f.write(f"  ELEMENT_FILE   elem.g.dat\n")
        f.write(f"  WEIGHT_FILE    weights.g.dat\n")
        f.write(f"  INDEX_FILE     index.g.dat\n")
        f.write("\n")
        
        f.write("LOCAL_MESH\n")
        f.write(f"  NODE_FILE      node.l.dat\n")
        f.write(f"  ELEMENT_FILE   elem.l.dat\n")
        f.write(f"  WEIGHT_FILE    weights.l.dat\n")
        f.write("\n")
        
        # Virtual mesh (if applicable)
        f.write("VIRTUAL_MESH\n")
        f.write(f"  NODE_FILE      node.v.dat\n")
        f.write(f"  ELEMENT_FILE   elem.v.dat\n")
        f.write("\n")
        
        # Boundary conditions
        f.write("# Boundary conditions\n")
        f.write("BOUNDARY_CONDITIONS\n")
        f.write(f"  GLOBAL_BC_FILE   bc.g.dat\n")
        f.write(f"  LOCAL_BC_FILE    bc.l.dat\n")
        f.write(f"  LOAD_FILE        load.dat\n")
        f.write("\n")
        
        # Material properties
        f.write("# Material properties\n")
        f.write("MATERIAL\n")
        f.write(f"  YOUNGS_MODULUS   {mp.EE:.15e}\n")
        f.write(f"  POISSON_RATIO    {mp.Nu:.15e}\n")
        f.write(f"  DENSITY          {mp.Rho:.15e}\n")
        f.write(f"  YIELD_STRESS     {mp.SigmaY0:.15e}\n")
        f.write("\n")
        
        # Time integration parameters (HHT method)
        f.write("# Time integration (HHT method)\n")
        f.write("TIME_INTEGRATION\n")
        f.write(f"  ALPHA   {sp.Alpha:.15e}\n")
        f.write(f"  BETA    {sp.Beta:.15e}\n")
        f.write(f"  GAMMA   {sp.Gamma:.15e}\n")
        f.write("\n")
        
        # Rayleigh damping
        f.write("# Rayleigh damping\n")
        f.write("DAMPING\n")
        f.write(f"  ALPHA_L   {sp.Alpha_l:.15e}\n")
        f.write(f"  BETA_L    {sp.Beta_l:.15e}\n")
        f.write("\n")
        
        # Analysis parameters
        f.write("# Analysis parameters\n")
        f.write("ANALYSIS\n")
        f.write(f"  STEP           {step}\n")
        f.write(f"  RESTART        {REstart}\n")
        f.write(f"  INCREMENT      {sp.inc}\n")
        f.write(f"  NREFINE        {sp.nrefLlist}\n")
        f.write(f"  NGP            {sp.ngp}\n")
        f.write("\n")
        
        # Crack geometry
        f.write("# Crack geometry\n")
        f.write("CRACK\n")
        f.write(f"  RADIUS         {sp.c:.15e}\n")
        f.write(f"  THICKNESS      {sp.thi:.15e}\n")
        f.write(f"  VELOCITY       {sp.V:.15e}\n")
        f.write("\n")
        
        # Solver settings
        f.write("# Solver settings\n")
        f.write("SOLVER\n")
        f.write(f"  OPENMP_THREADS   {sp.OPENMP}\n")
        f.write(f"  ABORT            {sp.ABO}\n")
        f.write("\n")
        
        # Output settings
        f.write("# Output settings\n")
        f.write("OUTPUT\n")
        f.write(f"  FOLDER           {sp.OUTPUT_FOLDER}\n")
        f.write(f"  CALC_STEP        {sp.CALC_STEP}\n")
        f.write("\n")
        
        # Restart file (if applicable)
        if REstart == 1 and step > 0:
            prev_step = step - 1
            prev_str = step2str(prev_step)
            f.write("# Restart data\n")
            f.write("RESTART_FILES\n")
            f.write(f"  DISPLACEMENT_FILE   ../step{prev_str}/delta_u.dat\n")
            f.write(f"  VELOCITY_FILE       ../step{prev_str}/velocity.dat\n")
            f.write(f"  ACCELERATION_FILE   ../step{prev_str}/acceleration.dat\n")
            f.write("\n")
    
    logger.info("input.dat generated successfully")


def generate_virtual_mesh():
    """
    Generate virtual mesh for visualization or special purposes
    This is a placeholder - implement based on specific needs
    """
    logger.info("Generating virtual mesh (placeholder)")
    
    # Create simple virtual mesh
    # This could be used for visualization or output purposes
    nodes_v = np.array([[1, 0.0, 0.0, 0.0]])
    elems_v = np.array([[1, 1]])
    
    np.savetxt('node.v.dat', nodes_v, fmt=['%d', '%.15e', '%.15e', '%.15e'])
    np.savetxt('elem.v.dat', elems_v, fmt='%d')
    
    logger.info("Virtual mesh files generated")
