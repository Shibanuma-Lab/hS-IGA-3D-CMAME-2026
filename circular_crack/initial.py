import numpy as np
import os
from const import simulation_params as sp
from utils.logger import logger
from utils.step2str import step2str


def initial(step, local_mesh, global_mesh):
    """
    Initialize displacement, velocity, and acceleration fields for restart
    Reads previous step results and prepares initial conditions
    """
    logger.info(f"Initializing fields for restart at step {step}")
    
    if step == 0:
        logger.info("Step 0: No restart data needed")
        return None
    
    # Get previous step
    prev_step = step - 1
    prev_str = step2str(prev_step)
    
    # Path to previous results
    prev_path = f"../step{prev_str}"
    
    # Initialize data structures
    init_data = {
        'displacement': None,
        'velocity': None,
        'acceleration': None
    }
    
    # Read previous displacement
    disp_file = os.path.join(prev_path, 'delta_u.dat')
    if os.path.exists(disp_file):
        init_data['displacement'] = np.loadtxt(disp_file)
        logger.info(f"Loaded displacement from {disp_file}")
    else:
        logger.warning(f"Displacement file not found: {disp_file}")
        # Initialize with zeros
        n_dof = _get_total_dof(global_mesh, local_mesh)
        init_data['displacement'] = np.zeros(n_dof)
    
    # Read previous velocity
    vel_file = os.path.join(prev_path, 'velocity.dat')
    if os.path.exists(vel_file):
        init_data['velocity'] = np.loadtxt(vel_file)
        logger.info(f"Loaded velocity from {vel_file}")
    else:
        logger.warning(f"Velocity file not found: {vel_file}")
        n_dof = _get_total_dof(global_mesh, local_mesh)
        init_data['velocity'] = np.zeros(n_dof)
    
    # Read previous acceleration
    acc_file = os.path.join(prev_path, 'acceleration.dat')
    if os.path.exists(acc_file):
        init_data['acceleration'] = np.loadtxt(acc_file)
        logger.info(f"Loaded acceleration from {acc_file}")
    else:
        logger.warning(f"Acceleration file not found: {acc_file}")
        n_dof = _get_total_dof(global_mesh, local_mesh)
        init_data['acceleration'] = np.zeros(n_dof)
    
    # Write initial condition files for current step
    _write_initial_files(init_data)
    
    logger.info("Initial conditions prepared successfully")
    return init_data


def _get_total_dof(global_mesh, local_mesh):
    """
    Calculate total degrees of freedom
    """
    n_nodes_g = len(global_mesh.node_g) if global_mesh.node_g is not None else 0
    n_nodes_l = len(local_mesh.node_l) if local_mesh.node_l is not None else 0
    
    # 3 DOF per node (ux, uy, uz)
    return 3 * (n_nodes_g + n_nodes_l)


def _write_initial_files(init_data):
    """
    Write initial condition files
    """
    if init_data['displacement'] is not None:
        np.savetxt('delta_u_init.dat', init_data['displacement'], fmt='%.15e')
    
    if init_data['velocity'] is not None:
        np.savetxt('velocity_init.dat', init_data['velocity'], fmt='%.15e')
    
    if init_data['acceleration'] is not None:
        np.savetxt('acceleration_init.dat', init_data['acceleration'], fmt='%.15e')


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
