"""
FEM data loader for reading displacement, velocity, and acceleration from .mat files
"""

import numpy as np
import os
from scipy.io import loadmat
from utils.logger import logger


class FEMDataLoader:
    """
    Load FEM data from .mat files for boundary condition interpolation
    """
    
    def __init__(self, fem_data_folder=None):
        """
        Initialize FEM data loader
        
        Args:
            fem_data_folder: Path to folder containing .mat files
        """
        if fem_data_folder is None:
            # Use path relative to this file's location
            current_dir = os.path.dirname(os.path.abspath(__file__))
            fem_data_folder = os.path.join(current_dir, 'data', 'FEMdata')
        self.fem_data_folder = fem_data_folder
        self.data = None
        self.elemGFEM = None
        self.nodeGFEM = None
        self.disG2DAllMa2D = None
        self.velG2DAllMa2D = None
        self.acceG2DAllMa2D = None
        
    def load_data(self, velocity):
        """
        Load FEM data for given crack velocity
        
        Args:
            velocity: Crack propagation velocity (e.g., 1000)
        
        Returns:
            True if successful, False otherwise
        """
        # Convert to int to avoid float formatting issues
        velocity = int(velocity)
        mat_file = os.path.join(self.fem_data_folder, f'FEM_v_{velocity}_uva.mat')
        
        if not os.path.exists(mat_file):
            logger.warning(f"FEM data file not found: {mat_file}")
            return False
        
        try:
            logger.info(f"Loading FEM data from {mat_file}")
            mat_data = loadmat(mat_file)
            
            # Extract data from Expression1 structure
            expr1 = mat_data['Expression1'][0, 0]
            self.elemGFEM = expr1['elemGFEM']
            self.nodeGFEM = expr1['nodeGFEM']
            self.disG2DAllMa2D = expr1['disG2DAllMa2D']
            self.velG2DAllMa2D = expr1['velG2DAllMa2D'] if 'velG2DAllMa2D' in expr1.dtype.names else None
            self.acceG2DAllMa2D = expr1['acceG2DAllMa2D'] if 'acceG2DAllMa2D' in expr1.dtype.names else None
            
            logger.info(f"FEM data loaded successfully:")
            logger.info(f"  Elements: {self.elemGFEM.shape}")
            logger.info(f"  Nodes: {self.nodeGFEM.shape}")
            logger.info(f"  Displacement: {self.disG2DAllMa2D.shape}")
            logger.info(f"  Velocity: {self.velG2DAllMa2D.shape}")
            logger.info(f"  Acceleration: {self.acceG2DAllMa2D.shape}")
            
            return True
            
        except Exception as e:
            logger.error(f"Error loading FEM data: {e}")
            return False
    
    def get_displacement_at_step(self, step):
        """
        Get displacement field at given step
        
        Args:
            step: Time step number (0-based)
        
        Returns:
            Displacement array of shape (n_nodes, 2) or None if not available
        """
        if self.disG2DAllMa2D is None:
            return None
        
        if step >= self.disG2DAllMa2D.shape[0]:
            logger.warning(f"Step {step} exceeds available FEM data ({self.disG2DAllMa2D.shape[0]} steps)")
            return None
        
        return self.disG2DAllMa2D[step]
    
    def get_velocity_at_step(self, step):
        """
        Get velocity field at given step
        
        Args:
            step: Time step number (0-based)
        
        Returns:
            Velocity array of shape (n_nodes, 2) or None if not available
        """
        if self.velG2DAllMa2D is None:
            return None
        
        if step >= self.velG2DAllMa2D.shape[0]:
            logger.warning(f"Step {step} exceeds available FEM data ({self.velG2DAllMa2D.shape[0]} steps)")
            return None
        
        return self.velG2DAllMa2D[step]
    
    def get_acceleration_at_step(self, step):
        """
        Get acceleration field at given step
        
        Args:
            step: Time step number (0-based)
        
        Returns:
            Acceleration array of shape (n_nodes, 2) or None if not available
        """
        if self.acceG2DAllMa2D is None:
            return None
        
        if step >= self.acceG2DAllMa2D.shape[0]:
            logger.warning(f"Step {step} exceeds available FEM data ({self.acceG2DAllMa2D.shape[0]} steps)")
            return None
        
        return self.acceG2DAllMa2D[step]
    
    def interpolate_2d_to_3d(self, values_2d, nodeG_3d):
        """
        Interpolate 2D FEM values to 3D IGA nodes with cylindrical coordinate transformation
        
        FEM data is in cylindrical coordinates (r, z):
        - nodeGFEM[:, 0] = r (radial position)
        - nodeGFEM[:, 1] = z (axial position)
        - values_2d[:, 0] = u_r (radial displacement)
        - values_2d[:, 1] = u_z (axial displacement)
        
        Transformation to 3D Cartesian:
        - For each 3D node at (x, y, z):
          1. Convert to cylindrical: r = sqrt(x^2 + y^2), theta = atan2(y, x)
          2. Interpolate (u_r, u_z) at (r, z) from FEM data
          3. Convert to Cartesian: u_x = u_r * cos(theta), u_y = u_r * sin(theta), u_z = u_z
        
        Args:
            values_2d: (n_fem_nodes, 2) array - [u_r, u_z] in cylindrical coordinates
            nodeG_3d: (n_iga_nodes, 3) array - [x, y, z] in Cartesian coordinates
        
        Returns:
            (n_iga_nodes, 3) array - [u_x, u_y, u_z] in Cartesian coordinates
        """
        from scipy.interpolate import LinearNDInterpolator
        
        # Create interpolators for radial and axial components
        # FEM nodes are in (r, z) coordinates
        interp_ur = LinearNDInterpolator(self.nodeGFEM, values_2d[:, 0])  # u_r
        interp_uz = LinearNDInterpolator(self.nodeGFEM, values_2d[:, 1])  # u_z
        
        # Convert 3D IGA nodes to cylindrical coordinates
        x = nodeG_3d[:, 0]
        y = nodeG_3d[:, 1]
        z = nodeG_3d[:, 2]
        
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        
        # Create (r, z) points for interpolation
        rz_points = np.column_stack([r, z])
        
        # Interpolate u_r and u_z
        u_r = interp_ur(rz_points)
        u_z = interp_uz(rz_points)
        
        # Replace NaN with 0 (points outside convex hull)
        u_r = np.nan_to_num(u_r)
        u_z = np.nan_to_num(u_z)
        
        # Convert cylindrical displacements to Cartesian coordinates
        values_3d = np.zeros((len(nodeG_3d), 3))
        values_3d[:, 0] = u_r * np.cos(theta)  # u_x = u_r * cos(theta)
        values_3d[:, 1] = u_r * np.sin(theta)  # u_y = u_r * sin(theta)
        values_3d[:, 2] = u_z                   # u_z = u_z
        
        return values_3d
