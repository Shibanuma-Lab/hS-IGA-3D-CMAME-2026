import numpy as np
import os
from const import simulation_params as sp, const_local_mesh as clm
from utils.logger import logger
from utils.step2str import step2str


def initial(step, local_mesh, global_mesh):
    """
    Initialize displacement, velocity, and acceleration fields for restart
    Based on Mathematica's initial[step] function
    
    Generates 6 files:
    - init.u.g.dat, init.v.g.dat, init.a.g.dat (global mesh)
    - init.u.l.dat, init.v.l.dat, init.a.l.dat (local mesh)
    """
    logger.info(f"Initializing fields for restart at step {step}")
    
    if step == 0:
        logger.info("Step 0: No initialization needed")
        return None
    
    # Get previous step
    prev_step = step - 1
    prev_str = step2str(prev_step)
    
    # Get project root directory and construct path to results
    current_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(current_dir, "results")
    prev_path = os.path.join(results_dir, f"step{prev_str}", "log")
    
    # Read previous step results from log directory
    disG, velG, acceG, disL, velL, acceL = _read_previous_results(prev_path, global_mesh, local_mesh)
    
    # Process initial conditions based on step
    if 1 <= step <= clm.aL:
        # For steps 1 to aL: use previous results as-is
        disiniG = disG
        veliniG = velG
        acceiniG = acceG
        disiniL = disL
        veliniL = velL
        acceiniL = acceL
        logger.info(f"Step {step} <= aL: Using previous results directly")
    
    elif step > clm.aL:
        # For steps > aL: global mesh unchanged, but local mesh nodes at r=rmax set to zero
        disiniG = disG
        veliniG = velG
        acceiniG = acceG
        
        # For local mesh: nodes where Mod[#, nLr+1] == 0 (r = rmax) set to zero
        nLr = clm.aL + clm.lL
        nnmL = len(local_mesh.nodeL)
        
        disiniL = np.array([disL[i] if (i+1) % (nLr + 1) != 0 else np.array([0., 0., 0.]) 
                           for i in range(nnmL)])
        veliniL = np.array([velL[i] if (i+1) % (nLr + 1) != 0 else np.array([0., 0., 0.]) 
                           for i in range(nnmL)])
        acceiniL = np.array([acceL[i] if (i+1) % (nLr + 1) != 0 else np.array([0., 0., 0.]) 
                            for i in range(nnmL)])
        logger.info(f"Step {step} > aL: Reset r=rmax nodes to zero in local mesh")
    
    # Write initial condition files
    _write_initial_files(disiniG, veliniG, acceiniG, disiniL, veliniL, acceiniL)
    
    logger.info("Initial condition files written successfully")
    return {
        'disG': disiniG,
        'velG': veliniG,
        'acceG': acceiniG,
        'disL': disiniL,
        'velL': veliniL,
        'acceL': acceiniL
    }


def _read_previous_results(prev_path, global_mesh, local_mesh):
    """
    Read displacement, velocity, acceleration from previous step's log directory
    Format: First line is node count, subsequent lines: nodeID x y z
    """
    nnmG = len(global_mesh.nodeG)
    nnmL = len(local_mesh.nodeL)
    
    # Try to read from log directory
    try:
        # Global mesh results
        disG_file = os.path.join(prev_path, 'u.g.dat')
        velG_file = os.path.join(prev_path, 'v.g.dat')
        acceG_file = os.path.join(prev_path, 'a.g.dat')
        
        if os.path.exists(disG_file):
            disG_data = np.loadtxt(disG_file, skiprows=1)
            disG = disG_data[:, 1:4]  # Skip node ID, take columns 2-4
            logger.info(f"Loaded global displacement: {disG.shape}")
        else:
            logger.warning(f"Global displacement file not found, using zeros")
            disG = np.zeros((nnmG, 3))
        
        if os.path.exists(velG_file):
            velG_data = np.loadtxt(velG_file, skiprows=1)
            velG = velG_data[:, 1:4]
            logger.info(f"Loaded global velocity: {velG.shape}")
        else:
            logger.warning(f"Global velocity file not found, using zeros")
            velG = np.zeros((nnmG, 3))
        
        if os.path.exists(acceG_file):
            acceG_data = np.loadtxt(acceG_file, skiprows=1)
            acceG = acceG_data[:, 1:4]
            logger.info(f"Loaded global acceleration: {acceG.shape}")
        else:
            logger.warning(f"Global acceleration file not found, using zeros")
            acceG = np.zeros((nnmG, 3))
        
        # Local mesh results
        disL_file = os.path.join(prev_path, 'u.l.dat')
        velL_file = os.path.join(prev_path, 'v.l.dat')
        acceL_file = os.path.join(prev_path, 'a.l.dat')
        
        if os.path.exists(disL_file):
            disL_data = np.loadtxt(disL_file, skiprows=1)
            disL = disL_data[:, 1:4]
            logger.info(f"Loaded local displacement: {disL.shape}")
        else:
            logger.warning(f"Local displacement file not found, using zeros")
            disL = np.zeros((nnmL, 3))
        
        if os.path.exists(velL_file):
            velL_data = np.loadtxt(velL_file, skiprows=1)
            velL = velL_data[:, 1:4]
            logger.info(f"Loaded local velocity: {velL.shape}")
        else:
            logger.warning(f"Local velocity file not found, using zeros")
            velL = np.zeros((nnmL, 3))
        
        if os.path.exists(acceL_file):
            acceL_data = np.loadtxt(acceL_file, skiprows=1)
            acceL = acceL_data[:, 1:4]
            logger.info(f"Loaded local acceleration: {acceL.shape}")
        else:
            logger.warning(f"Local acceleration file not found, using zeros")
            acceL = np.zeros((nnmL, 3))
    
    except Exception as e:
        logger.error(f"Error reading previous results: {e}")
        logger.info("Initializing with zeros")
        disG = np.zeros((nnmG, 3))
        velG = np.zeros((nnmG, 3))
        acceG = np.zeros((nnmG, 3))
        disL = np.zeros((nnmL, 3))
        velL = np.zeros((nnmL, 3))
        acceL = np.zeros((nnmL, 3))
    
    return disG, velG, acceG, disL, velL, acceL


def _write_initial_files(disG, velG, acceG, disL, velL, acceL):
    """
    Write initial condition files in Mathematica format
    Format: Line 1: node count
            Line 2+: nodeID \t x \t y \t z (15 significant digits)
    """
    
    from utils.format_output import format_real
    
    def write_node_data(filename, data):
        """Write node data in format: nodeID x y z"""
        lines = [str(len(data))]
        for i, (x, y, z) in enumerate(data, start=1):
            line = f"{i}\t{format_real(x)}\t{format_real(y)}\t{format_real(z)}"
            lines.append(line)
        with open(filename, 'w') as f:
            f.write('\n'.join(lines))
    
    # Write global mesh initial conditions
    write_node_data('init.u.g.dat', disG)
    write_node_data('init.v.g.dat', velG)
    write_node_data('init.a.g.dat', acceG)
    logger.info("Global initial condition files written")
    
    # Write local mesh initial conditions
    write_node_data('init.u.l.dat', disL)
    write_node_data('init.v.l.dat', velL)
    write_node_data('init.a.l.dat', acceL)
    logger.info("Local initial condition files written")


def generate_zero_initial_conditions(global_mesh, local_mesh):
    """
    Generate zero initial conditions for the first step
    """
    logger.info("Generating zero initial conditions")
    
    n_dof = _get_total_dof(global_mesh, local_mesh)
    
    init_data = {
        'displacement': np.zeros(n_dof),
        'velocity': np.zeros(n_dof),
        'acceleration': np.zeros(n_dof)
    }
    
    _write_initial_files(init_data)
    
    return init_data


def apply_initial_velocity(init_data, node_ids, velocity_vector):
    """
    Apply initial velocity to specific nodes
    Useful for dynamic problems with prescribed initial conditions
    
    Parameters:
    -----------
    init_data : dict
        Dictionary containing displacement, velocity, acceleration
    node_ids : array-like
        Node IDs to apply velocity
    velocity_vector : array-like
        Velocity vector [vx, vy, vz] for each node
    """
    for i, node_id in enumerate(node_ids):
        # Convert node ID to DOF indices
        dof_start = 3 * (node_id - 1)
        init_data['velocity'][dof_start:dof_start+3] = velocity_vector[i]
    
    return init_data


def calculate_crack_growth_rate(step):
    """
    Calculate crack growth rate at current step
    Used for dynamic crack propagation
    """
    if step == 0:
        return 0.0
    
    # Crack length at current and previous step
    a_current = sp.c * step / sp.stepall
    a_previous = sp.c * (step - 1) / sp.stepall
    
    # Time increment (calculated from crack velocity)
    dt = (a_current - a_previous) / sp.V
    
    return sp.V, dt
