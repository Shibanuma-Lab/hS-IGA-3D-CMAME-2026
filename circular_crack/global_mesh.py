import numpy as np
from const import const_global_mesh as cgm, const_local_mesh as clm, simulation_params as sp
from utils.logger import logger


class GlobalMesh:
    """
    Global mesh generator for circular crack problem
    Based on Mathematica implementation
    """
    
    def __init__(self, step):
        self.step = step
        self.m0 = cgm.m0
        self.nx1 = cgm.nx1
        self.ny1 = cgm.ny1
        self.rGL = cgm.rGL
        
        # Local mesh parameters
        self.hL = clm.hL
        self.aL = clm.aL
        self.lL = clm.lL
        self.HL = clm.HL
        
        # Simulation parameters
        self.c = sp.c  # crack radius
        self.thi = sp.thi  # thickness
        
        # Mesh data
        self.node_g = None
        self.elem_g = None
        self.weights_g = None
        self.index_g = None
        
        # B-spline degree
        self.deg = 3
        
    def make_global_mesh(self):
        """
        Generate global mesh with graded elements
        """
        logger.info("Generating global mesh nodes and elements")
        
        # Calculate global element size at crack surface
        hG = self.rGL * self.hL
        
        # Calculate crack length at current step
        a = self.c * self.step / sp.stepall
        
        # Generate mesh sizes with grading
        # m0 levels of refinement from hG to larger elements
        mesh_sizes = [hG * (2 ** i) for i in range(self.m0)]
        
        # X direction: crack region + additional region
        x_coords = self._generate_graded_coords_x(a, hG, mesh_sizes)
        
        # Y direction: through thickness
        y_coords = self._generate_graded_coords_y(mesh_sizes)
        
        # Z direction: perpendicular to crack plane
        z_coords = self._generate_graded_coords_z(mesh_sizes)
        
        # Generate nodes as tensor product
        nodes = []
        node_id = 1
        for k, z in enumerate(z_coords):
            for j, y in enumerate(y_coords):
                for i, x in enumerate(x_coords):
                    nodes.append([node_id, x, y, z])
                    node_id += 1
        
        self.node_g = np.array(nodes)
        
        # Generate elements (NURBS patches)
        self._generate_elements(len(x_coords), len(y_coords), len(z_coords))
        
        # Generate weights (all 1.0 for B-spline, can be modified for NURBS)
        self._generate_weights()
        
        # Generate control point indices
        self._generate_indices()
        
        logger.info(f"Global mesh generated: {len(self.node_g)} nodes, {len(self.elem_g)} elements")
        
    def _generate_graded_coords_x(self, a, hG, mesh_sizes):
        """
        Generate graded coordinates in X direction
        """
        coords = [0.0]
        
        # Fine mesh near crack tip
        x = 0.0
        # Region with finest mesh (crack length)
        n_fine = int(a / hG)
        for i in range(n_fine):
            x += hG
            coords.append(x)
        
        # Graded mesh (coarser elements away from crack)
        for level, h in enumerate(mesh_sizes[1:], start=1):
            n_elements = self.nx1 // (2 ** level)
            for i in range(max(1, n_elements)):
                x += h
                coords.append(x)
        
        return np.array(coords)
    
    def _generate_graded_coords_y(self, mesh_sizes):
        """
        Generate graded coordinates in Y direction (thickness)
        """
        coords = [-self.thi / 2]
        
        # Symmetric grading from center to edges
        y = -self.thi / 2
        
        # Fine elements in the middle
        h = mesh_sizes[0]
        n_fine = self.ny1
        
        for i in range(n_fine):
            y += h
            coords.append(y)
        
        # Coarser elements toward edges
        for h in mesh_sizes[1:]:
            for i in range(max(1, self.ny1 // 2)):
                y += h
                if y <= self.thi / 2:
                    coords.append(y)
        
        # Ensure we reach the upper surface
        if coords[-1] < self.thi / 2:
            coords.append(self.thi / 2)
        
        return np.array(coords)
    
    def _generate_graded_coords_z(self, mesh_sizes):
        """
        Generate graded coordinates in Z direction
        """
        coords = [0.0]
        
        z = 0.0
        # Start with fine mesh
        h = mesh_sizes[0]
        
        # Create graded mesh in Z direction
        for level, h in enumerate(mesh_sizes):
            n_elements = max(1, self.ny1 // (2 ** level))
            for i in range(n_elements):
                z += h
                coords.append(z)
        
        return np.array(coords)
    
    def _generate_elements(self, nx, ny, nz):
        """
        Generate element connectivity for NURBS/B-spline patches
        For degree p, each element uses (p+1)^3 control points
        """
        elements = []
        elem_id = 1
        p = self.deg  # degree
        
        # For B-spline, elements span between knot spans
        for k in range(nz - p):
            for j in range(ny - p):
                for i in range(nx - p):
                    # Control points for this element
                    ctrl_pts = []
                    for kk in range(p + 1):
                        for jj in range(p + 1):
                            for ii in range(p + 1):
                                node_idx = (k + kk) * nx * ny + (j + jj) * nx + (i + ii) + 1
                                ctrl_pts.append(node_idx)
                    
                    elements.append([elem_id] + ctrl_pts)
                    elem_id += 1
        
        self.elem_g = np.array(elements, dtype=int)
    
    def _generate_weights(self):
        """
        Generate weights for NURBS (all 1.0 for B-spline)
        """
        n_nodes = len(self.node_g)
        self.weights_g = np.ones((n_nodes, 2))
        self.weights_g[:, 0] = np.arange(1, n_nodes + 1)  # Node ID
        # weights_g[:, 1] are all 1.0
    
    def _generate_indices(self):
        """
        Generate control point indices
        """
        n_nodes = len(self.node_g)
        self.index_g = np.arange(1, n_nodes + 1).reshape(-1, 1)
        self.index_g = np.column_stack([self.index_g, self.index_g])
    
    def generate(self):
        """
        Write mesh data to files
        """
        logger.info("Writing global mesh files")
        
        # Write node file
        np.savetxt('node.g.dat', self.node_g, fmt=['%d', '%.15e', '%.15e', '%.15e'])
        
        # Write element file
        fmt_elem = ['%d'] + ['%d'] * (self.elem_g.shape[1] - 1)
        np.savetxt('elem.g.dat', self.elem_g, fmt=fmt_elem)
        
        # Write weights file
        np.savetxt('weights.g.dat', self.weights_g, fmt=['%d', '%.15e'])
        
        # Write index file
        np.savetxt('index.g.dat', self.index_g, fmt=['%d', '%d'])
        
        logger.info("Global mesh files written successfully")
