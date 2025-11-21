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
        Interpolate 2D FEM values to 3D IGA nodes using exact Mathematica algorithm
        
        This implements the exact same interpolation as Mathematica's boundaryebcG:
        1. For each IGA node, find containing FEM element using bounding box search
        2. Solve for parametric coordinates (ξ, η) within element
        3. Use bilinear shape functions N1, N2, N3, N4 to interpolate displacement
        4. Transform from cylindrical (u_r, u_z) to Cartesian (u_x, u_y, u_z)
        
        FEM data is in cylindrical coordinates (r, z):
        - nodeGFEM[:, 0] = r (radial position)
        - nodeGFEM[:, 1] = z (axial position)
        - values_2d[:, 0] = u_r (radial displacement)
        - values_2d[:, 1] = u_z (axial displacement)
        
        Args:
            values_2d: (n_fem_nodes, 2) array - [u_r, u_z] in cylindrical coordinates
            nodeG_3d: (n_iga_nodes, 3) array - [x, y, z] in Cartesian coordinates
        
        Returns:
            (n_iga_nodes, 3) array - [u_x, u_y, u_z] in Cartesian coordinates
        """
        from scipy.optimize import fsolve
        
        # Bilinear shape functions (matches Mathematica exactly)
        def N1(xi, eta):
            return 0.25 * (1.0 - xi) * (1.0 - eta)
        
        def N2(xi, eta):
            return 0.25 * (1.0 + xi) * (1.0 - eta)
        
        def N3(xi, eta):
            return 0.25 * (1.0 + xi) * (1.0 + eta)
        
        def N4(xi, eta):
            return 0.25 * (1.0 - xi) * (1.0 + eta)
        
        # Pre-compute element bounding boxes for faster search
        neGa = len(self.elemGFEM)
        elGaposGa = []
        rmaxGa = np.zeros(neGa)
        rminGa = np.zeros(neGa)
        zmaxGa = np.zeros(neGa)
        zminGa = np.zeros(neGa)
        
        for eG in range(neGa):
            # Get element node positions (convert from 1-based to 0-based indexing)
            elem_node_indices = self.elemGFEM[eG] - 1
            elem_nodes = self.nodeGFEM[elem_node_indices]
            elGaposGa.append(elem_nodes)
            rmaxGa[eG] = np.max(elem_nodes[:, 0])
            rminGa[eG] = np.min(elem_nodes[:, 0])
            zmaxGa[eG] = np.max(elem_nodes[:, 1])
            zminGa[eG] = np.min(elem_nodes[:, 1])
        
        def solve_xi_eta_Ga(eGa, rz_point):
            """
            Solve for parametric coordinates (ξ, η) given (r, z) point in element
            Matches Mathematica's solve\[Xi]\[Eta]Ga function
            """
            xy = elGaposGa[eGa]
            r_target, z_target = rz_point
            
            # Check if element is aligned rectangular (optimization from Mathematica)
            if (np.abs(xy[0, 0] - xy[3, 0]) < 1e-10 and 
                np.abs(xy[1, 0] - xy[2, 0]) < 1e-10 and
                np.abs(xy[0, 1] - xy[1, 1]) < 1e-10 and
                np.abs(xy[2, 1] - xy[3, 1]) < 1e-10):
                # Rectangular element aligned with axes
                r_center = 0.5 * (xy[2, 0] + xy[0, 0])
                z_center = 0.5 * (xy[2, 1] + xy[0, 1])
                r_half = 0.5 * (xy[2, 0] - xy[0, 0])
                z_half = 0.5 * (xy[2, 1] - xy[0, 1])
                
                # Avoid division by zero
                if abs(r_half) < 1e-15:
                    xi = 0.0
                else:
                    xi = (r_target - r_center) / r_half
                
                if abs(z_half) < 1e-15:
                    eta = 0.0
                else:
                    eta = (z_target - z_center) / z_half
                
                return xi, eta
            else:
                # General case: solve nonlinear system using Newton-Raphson
                def residual(params):
                    xi, eta = params
                    r_interp = (N1(xi, eta) * xy[0, 0] + 
                               N2(xi, eta) * xy[1, 0] + 
                               N3(xi, eta) * xy[2, 0] + 
                               N4(xi, eta) * xy[3, 0])
                    z_interp = (N1(xi, eta) * xy[0, 1] + 
                               N2(xi, eta) * xy[1, 1] + 
                               N3(xi, eta) * xy[2, 1] + 
                               N4(xi, eta) * xy[3, 1])
                    return [r_interp - r_target, z_interp - z_target]
                
                # Try multiple initial guesses if first one fails
                initial_guesses = [
                    [0.0, 0.0],      # Center
                    [-0.5, -0.5],    # Lower-left quadrant
                    [0.5, -0.5],     # Lower-right quadrant
                    [-0.5, 0.5],     # Upper-left quadrant
                    [0.5, 0.5],      # Upper-right quadrant
                ]
                
                for guess in initial_guesses:
                    try:
                        solution = fsolve(residual, guess, full_output=True)
                        xi, eta = solution[0]
                        info = solution[1]
                        
                        # Check if solution converged and is valid
                        if info['fvec'][0]**2 + info['fvec'][1]**2 < 1e-10:
                            return xi, eta
                    except:
                        continue
                
                # If all guesses fail, return center (0, 0)
                return 0.0, 0.0
        
        def get_xi_eta_GaeG(node_xyz):
            """
            Find containing element and parametric coordinates for 3D node
            Matches Mathematica's get\[Xi]\[Eta]GaeG function
            """
            r = np.sqrt(node_xyz[0]**2 + node_xyz[1]**2)
            z = node_xyz[2]
            
            # Find candidate elements using bounding box (matches Mathematica's eGaabb)
            eGaabb = []
            for eG in range(neGa):
                if (rmaxGa[eG] >= r and rminGa[eG] <= r and
                    zmaxGa[eG] >= z and zminGa[eG] <= z):
                    eGaabb.append(eG)
            
            if len(eGaabb) == 0:
                return None, None
            
            # Try each candidate element
            best_eGa = None
            best_xi_eta = None
            best_dist = float('inf')
            
            for eGa in eGaabb:
                xi, eta = solve_xi_eta_Ga(eGa, [r, z])
                
                # Calculate distance from valid range [-1, 1]
                dist = 0.0
                if xi < -1.0:
                    dist += (-1.0 - xi)**2
                elif xi > 1.0:
                    dist += (xi - 1.0)**2
                if eta < -1.0:
                    dist += (-1.0 - eta)**2
                elif eta > 1.0:
                    dist += (eta - 1.0)**2
                
                # Accept if within tolerance (allow small overshoot for boundary nodes)
                tolerance = 0.01  # Allow 1% overshoot for numerical errors
                if dist < tolerance**2:
                    if dist < best_dist:
                        best_eGa = eGa
                        best_xi_eta = (xi, eta)
                        best_dist = dist
                        
                        # Perfect match - return immediately
                        if dist == 0.0:
                            return eGa, (xi, eta)
            
            # Return best match even if slightly outside
            if best_eGa is not None:
                return best_eGa, best_xi_eta
            
            return None, None
        
        def disinter(eG, xi_eta, disp_data):
            """
            Interpolate displacement at parametric coordinates within element
            Matches Mathematica's disinter function exactly
            """
            xi, eta = xi_eta
            nn = [N1(xi, eta), N2(xi, eta), N3(xi, eta), N4(xi, eta)]
            
            # Get displacement at element nodes (convert from 1-based to 0-based)
            elem_node_indices = self.elemGFEM[eG] - 1
            disGn = disp_data[elem_node_indices]
            
            # Interpolate u_r and u_z
            u_r = sum(nn[i] * disGn[i, 0] for i in range(4))
            u_z = sum(nn[i] * disGn[i, 1] for i in range(4))
            
            return u_r, u_z
        
        # Interpolate for each IGA node
        values_3d = np.zeros((len(nodeG_3d), 3))
        
        for i, node_xyz in enumerate(nodeG_3d):
            # Find containing element and parametric coordinates
            eG, xi_eta = get_xi_eta_GaeG(node_xyz)
            
            if eG is None:
                # Point not found in any element (outside FEM domain)
                values_3d[i] = [0.0, 0.0, 0.0]
                continue
            
            # Interpolate displacement in cylindrical coordinates
            u_r, u_z = disinter(eG, xi_eta, values_2d)
            
            # Transform to Cartesian coordinates (matches Mathematica's getbc function)
            x, y = node_xyz[0], node_xyz[1]
            if x == 0.0:
                theta = np.pi / 2.0
            else:
                theta = np.arctan(y / x)
            
            values_3d[i, 0] = u_r * np.cos(theta)  # u_x
            values_3d[i, 1] = u_r * np.sin(theta)  # u_y
            values_3d[i, 2] = u_z                   # u_z
        
        return values_3d
