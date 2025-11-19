import numpy as np
from const import simulation_params as sp, material_property as mp
from utils.logger import logger


class Boundary:
    """
    Boundary condition handler for circular crack problem
    Defines displacement BCs, loads, and constraints
    """
    
    def __init__(self, local_mesh, global_mesh):
        self.local_mesh = local_mesh
        self.global_mesh = global_mesh
        
        # Boundary condition data
        self.bc_g = None  # Global mesh boundary conditions
        self.bc_l = None  # Local mesh boundary conditions
        self.load_data = None  # Applied loads
        
        # Material and loading parameters
        self.sigma_inf = mp.SigmaInfinity
        self.thi = sp.thi
        
    def define_boundary(self, local_mesh, global_mesh):
        """
        Define boundary conditions and loads
        Based on Mathematica boundary definition
        """
        logger.info("Defining boundary conditions and loads")
        
        # Define boundary conditions for global mesh
        self._define_global_bc(global_mesh)
        
        # Define boundary conditions for local mesh
        self._define_local_bc(local_mesh)
        
        # Define applied loads
        self._define_loads(global_mesh)
        
        logger.info(f"Boundary conditions defined: {len(self.bc_g)} global BCs, {len(self.bc_l)} local BCs")
    
    def _define_global_bc(self, global_mesh):
        """
        Define boundary conditions on global mesh
        Typically: fix nodes on back surface, apply symmetry conditions
        """
        bc_list = []
        
        nodes = global_mesh.node_g
        
        # Find nodes on boundaries
        # X = 0 plane: symmetry in X direction (ux = 0)
        x_min = np.min(nodes[:, 1])
        tol = 1e-10
        nodes_x_min = np.where(np.abs(nodes[:, 1] - x_min) < tol)[0]
        
        for node_idx in nodes_x_min:
            node_id = int(nodes[node_idx, 0])
            # Fix X displacement (symmetry)
            bc_list.append([node_id, 1, 0.0])  # DOF 1 (ux) = 0
        
        # Y = -thi/2 plane: symmetry in Y direction (uy = 0)
        y_min = np.min(nodes[:, 2])
        nodes_y_min = np.where(np.abs(nodes[:, 2] - y_min) < tol)[0]
        
        for node_idx in nodes_y_min:
            node_id = int(nodes[node_idx, 0])
            bc_list.append([node_id, 2, 0.0])  # DOF 2 (uy) = 0
        
        # Z = 0 plane: symmetry in Z direction (uz = 0)
        z_min = np.min(nodes[:, 3])
        nodes_z_min = np.where(np.abs(nodes[:, 3] - z_min) < tol)[0]
        
        for node_idx in nodes_z_min:
            node_id = int(nodes[node_idx, 0])
            bc_list.append([node_id, 3, 0.0])  # DOF 3 (uz) = 0
        
        self.bc_g = np.array(bc_list)
    
    def _define_local_bc(self, local_mesh):
        """
        Define boundary conditions on local mesh
        Local mesh needs to match global mesh at interface
        """
        bc_list = []
        
        nodes = local_mesh.node_l
        
        if nodes is None or len(nodes) == 0:
            self.bc_l = np.array([])
            return
        
        # Typically local mesh has constraints at the interface with global mesh
        # This depends on the coupling method (e.g., tied contact, constraint equations)
        
        # For now, apply symmetry conditions similar to global mesh
        tol = 1e-10
        
        # Y = -thi/2 plane: uy = 0
        y_min = np.min(nodes[:, 2])
        nodes_y_min = np.where(np.abs(nodes[:, 2] - y_min) < tol)[0]
        
        for node_idx in nodes_y_min:
            node_id = int(nodes[node_idx, 0])
            bc_list.append([node_id, 2, 0.0])
        
        # Z = 0 plane: uz = 0
        z_min = np.min(nodes[:, 3])
        nodes_z_min = np.where(np.abs(nodes[:, 3] - z_min) < tol)[0]
        
        for node_idx in nodes_z_min:
            node_id = int(nodes[node_idx, 0])
            bc_list.append([node_id, 3, 0.0])
        
        self.bc_l = np.array(bc_list) if bc_list else np.array([])
    
    def _define_loads(self, global_mesh):
        """
        Define applied loads (remote stress)
        Apply uniform tension stress at far field
        """
        load_list = []
        
        nodes = global_mesh.node_g
        
        # Find nodes on far field boundary (e.g., X = max)
        x_max = np.max(nodes[:, 1])
        tol = 1e-10
        nodes_x_max = np.where(np.abs(nodes[:, 1] - x_max) < tol)[0]
        
        # Calculate nodal force from uniform stress
        # For 3D problem: F = sigma * A / n_nodes
        # Area per node depends on mesh density
        
        # Simple approach: distribute total force equally
        if len(nodes_x_max) > 0:
            # This is a simplification - actual load distribution should consider
            # element shape functions and quadrature
            area_per_node = (self.thi * np.max(nodes[:, 3])) / len(nodes_x_max)
            force_per_node = self.sigma_inf * area_per_node
            
            for node_idx in nodes_x_max:
                node_id = int(nodes[node_idx, 0])
                # Apply force in X direction (DOF 1)
                load_list.append([node_id, 1, force_per_node])
        
        self.load_data = np.array(load_list) if load_list else np.array([])
    
    def generate(self):
        """
        Write boundary condition and load files
        """
        logger.info("Writing boundary condition files")
        
        # Write global BC file
        if self.bc_g is not None and len(self.bc_g) > 0:
            np.savetxt('bc.g.dat', self.bc_g, fmt=['%d', '%d', '%.15e'])
        
        # Write local BC file
        if self.bc_l is not None and len(self.bc_l) > 0:
            np.savetxt('bc.l.dat', self.bc_l, fmt=['%d', '%d', '%.15e'])
        
        # Write load file
        if self.load_data is not None and len(self.load_data) > 0:
            np.savetxt('load.dat', self.load_data, fmt=['%d', '%d', '%.15e'])
        
        logger.info("Boundary condition files written successfully")
