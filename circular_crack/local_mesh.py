import numpy as np
from const import const_local_mesh as clm, simulation_params as sp
from utils.logger import logger


class LocalMesh:
    """
    Local mesh generator for crack surface and ligament
    Creates refined mesh around crack tip
    """
    
    def __init__(self, step):
        self.step = step
        
        # Local mesh parameters
        self.hL = clm.hL  # Local element size
        self.aL = clm.aL  # Number of elements on crack surface
        self.lL = clm.lL  # Number of elements on ligament
        self.HL = clm.HL  # Number of elements in thickness
        
        # Simulation parameters
        self.c = sp.c  # Total crack radius
        self.thi = sp.thi  # Thickness
        
        # Mesh data
        self.node_l = None
        self.elem_l = None
        self.weights_l = None
        self.index_l = None
        
        # Crack geometry
        self.crack_length = None
        
        # B-spline degree
        self.deg = 3
        
    def make_local_mesh(self):
        """
        Generate local mesh for crack and ligament region
        """
        logger.info("Generating local mesh for crack region")
        
        # Calculate current crack length
        self.crack_length = self.c * self.step / sp.stepall
        
        # Generate coordinates for crack surface (polar-like mesh)
        r_coords, theta_coords, z_coords = self._generate_crack_coords()
        
        # Convert to Cartesian coordinates
        nodes = []
        node_id = 1
        
        for z in z_coords:
            for theta in theta_coords:
                for r in r_coords:
                    # Polar to Cartesian conversion
                    x = r * np.cos(theta)
                    y = r * np.sin(theta)
                    nodes.append([node_id, x, y, z])
                    node_id += 1
        
        self.node_l = np.array(nodes)
        
        # Generate elements
        nr = len(r_coords)
        nt = len(theta_coords)
        nz = len(z_coords)
        self._generate_elements(nr, nt, nz)
        
        # Generate weights
        self._generate_weights()
        
        # Generate indices
        self._generate_indices()
        
        logger.info(f"Local mesh generated: {len(self.node_l)} nodes, {len(self.elem_l)} elements")
    
    def _generate_crack_coords(self):
        """
        Generate coordinates for crack region in polar-like system
        """
        # Radial direction: from crack tip toward center
        # Fine mesh near crack tip
        r_coords = []
        r = 0.0
        for i in range(self.aL):
            r_coords.append(r)
            r += self.hL
        
        # Extend into ligament
        for i in range(self.lL):
            r_coords.append(r)
            r += self.hL
        
        r_coords = np.array(r_coords)
        
        # Angular direction: around crack front (0 to 2π for circular crack)
        # For circular crack, create angular mesh
        n_angular = max(8, int(2 * np.pi * self.crack_length / self.hL)) if self.crack_length > 0 else 8
        theta_coords = np.linspace(0, 2 * np.pi, n_angular + 1)[:-1]  # Remove duplicate at 2π
        
        # Through-thickness direction
        z_coords = np.linspace(-self.thi / 2, self.thi / 2, self.HL + 1)
        
        return r_coords, theta_coords, z_coords
    
    def _generate_elements(self, nr, nt, nz):
        """
        Generate element connectivity for local mesh
        """
        elements = []
        elem_id = 1
        p = self.deg
        
        # Handle periodic boundary in theta direction
        for k in range(nz - p):
            for j in range(nt):  # Periodic in theta
                for i in range(nr - p):
                    ctrl_pts = []
                    for kk in range(p + 1):
                        for jj in range(p + 1):
                            for ii in range(p + 1):
                                # Handle periodic wrapping in theta
                                j_idx = (j + jj) % nt
                                node_idx = (k + kk) * nr * nt + j_idx * nr + (i + ii) + 1
                                ctrl_pts.append(node_idx)
                    
                    elements.append([elem_id] + ctrl_pts)
                    elem_id += 1
        
        self.elem_l = np.array(elements, dtype=int)
    
    def _generate_weights(self):
        """
        Generate weights for NURBS
        For circular geometry, may need special weights
        """
        n_nodes = len(self.node_l)
        self.weights_l = np.ones((n_nodes, 2))
        self.weights_l[:, 0] = np.arange(1, n_nodes + 1)
        
        # For circular arc representation with NURBS, weights might need adjustment
        # Standard B-spline uses weight = 1.0
        # For exact circular arc, use weights like cos(theta/2)
        # Here we use simple B-spline approximation
    
    def _generate_indices(self):
        """
        Generate control point indices
        """
        n_nodes = len(self.node_l)
        self.index_l = np.arange(1, n_nodes + 1).reshape(-1, 1)
        self.index_l = np.column_stack([self.index_l, self.index_l])
    
    def generate(self):
        """
        Write local mesh data to files
        """
        logger.info("Writing local mesh files")
        
        # Write node file
        np.savetxt('node.l.dat', self.node_l, fmt=['%d', '%.15e', '%.15e', '%.15e'])
        
        # Write element file
        fmt_elem = ['%d'] + ['%d'] * (self.elem_l.shape[1] - 1)
        np.savetxt('elem.l.dat', self.elem_l, fmt=fmt_elem)
        
        # Write weights file
        np.savetxt('weights.l.dat', self.weights_l, fmt=['%d', '%.15e'])
        
        # Write index file (if needed by solver)
        # np.savetxt('index.l.dat', self.index_l, fmt=['%d', '%d'])
        
        logger.info("Local mesh files written successfully")
