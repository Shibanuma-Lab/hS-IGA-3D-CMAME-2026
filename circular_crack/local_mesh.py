"""
Local mesh generation for crack region
Direct translation from CircularCrackθL[] in function.txt
"""

import numpy as np
from const import const_local_mesh as clm, simulation_params as sp
from utils.logger import logger


class LocalMesh:
    """
    Local mesh generator for crack region
    Based on CircularCrackθL[step] from Mathematica
    """
    
    def __init__(self, step):
        self.step = step
        
        # Local mesh parameters
        self.hL = clm.hL
        self.aL = clm.aL
        self.lL = clm.lL
        self.HL = clm.HL
        self.d_theta = clm.d_theta
        
        # Mesh data
        self.nodeL = None
        self.elemL = None
        
        # Mesh dimensions
        self.nLr = None
        self.nLtheta = None
        self.nnmL = None
        self.nemL = None
        
    def make_local_mesh(self):
        """
        Generate local mesh
        Translation of CircularCrackθL[step]
        """
        logger.info(f"Generating local mesh for step {self.step}")
        
        # Helper function: set very small values to zero
        def ad(x):
            return 0.0 if abs(x) < 1e-9 * sp.WidthG else x
        
        # 1. Calculate mesh dimensions
        self.nLr = self.aL + self.lL
        self.nLtheta = round(90 / self.d_theta)
        
        # 2. Generate radial coordinates
        if self.step <= self.aL:
            # Initial crack: nodes from 0 to nLr
            nodeLr = self.hL * np.arange(self.nLr + 1)
        else:
            # Propagating crack: shift window
            nodeLr = self.hL * np.arange(self.step - self.aL, self.step + self.lL + 1)
        
        # 3. Generate through-thickness coordinates
        nodeLz = self.hL * np.arange(self.HL + 1)
        
        # 4. Generate angular coordinates (0 to 90 degrees)
        nodeLtheta = np.linspace(0, 0.5 * np.pi, self.nLtheta + 1)
        
        # 5. Create 3D nodes in cylindrical coordinates
        nodeL_list = []
        for z in nodeLz:
            for theta in nodeLtheta:
                for r in nodeLr:
                    x = ad(r * np.cos(theta))
                    y = ad(r * np.sin(theta))
                    z_val = ad(z)
                    nodeL_list.append([x, y, z_val])
        
        self.nodeL = np.array(nodeL_list)
        self.nnmL = len(self.nodeL)
        
        logger.info(f"Local nodes generated: {self.nnmL} nodes")
        logger.info(f"  nLr={self.nLr}, nLtheta={self.nLtheta}, HL={self.HL}")
        
        # 6. Generate element connectivity
        self._generate_elements()
        
        logger.info(f"Local mesh complete: {self.nemL} elements")
    
    def _generate_elements(self):
        """
        Generate 8-node hexahedral elements
        Translation from Mathematica eL calculation
        """
        NLrtheta = (self.nLr + 1) * (self.nLtheta + 1)
        
        # Initialize element connectivity (8 nodes per hex)
        eL = [[] for _ in range(8)]
        
        # Node numbering: varies fastest in r, then theta, then z
        # Element connectivity for 8-node hex:
        # Bottom face (z):     1-2-3-4
        # Top face (z+1):      5-6-7-8
        
        # eL[[1]]: node at (r, theta, z)
        for k in range(self.HL):
            for j in range(self.nLtheta):
                for i in range(self.nLr):
                    node_id = k * NLrtheta + j * (self.nLr + 1) + i + 1
                    eL[0].append(node_id)
        
        # eL[[2]]: node at (r+1, theta, z)
        for k in range(self.HL):
            for j in range(self.nLtheta):
                for i in range(1, self.nLr + 1):
                    node_id = k * NLrtheta + j * (self.nLr + 1) + i + 1
                    eL[1].append(node_id)
        
        # eL[[3]]: node at (r+1, theta+1, z)
        for k in range(self.HL):
            for j in range(1, self.nLtheta + 1):
                for i in range(1, self.nLr + 1):
                    node_id = k * NLrtheta + j * (self.nLr + 1) + i + 1
                    eL[2].append(node_id)
        
        # eL[[4]]: node at (r, theta+1, z)
        for k in range(self.HL):
            for j in range(1, self.nLtheta + 1):
                for i in range(self.nLr):
                    node_id = k * NLrtheta + j * (self.nLr + 1) + i + 1
                    eL[3].append(node_id)
        
        # eL[[5]]: node at (r, theta, z+1)
        for k in range(1, self.HL + 1):
            for j in range(self.nLtheta):
                for i in range(self.nLr):
                    node_id = k * NLrtheta + j * (self.nLr + 1) + i + 1
                    eL[4].append(node_id)
        
        # eL[[6]]: node at (r+1, theta, z+1)
        for k in range(1, self.HL + 1):
            for j in range(self.nLtheta):
                for i in range(1, self.nLr + 1):
                    node_id = k * NLrtheta + j * (self.nLr + 1) + i + 1
                    eL[5].append(node_id)
        
        # eL[[7]]: node at (r+1, theta+1, z+1)
        for k in range(1, self.HL + 1):
            for j in range(1, self.nLtheta + 1):
                for i in range(1, self.nLr + 1):
                    node_id = k * NLrtheta + j * (self.nLr + 1) + i + 1
                    eL[6].append(node_id)
        
        # eL[[8]]: node at (r, theta+1, z+1)
        for k in range(1, self.HL + 1):
            for j in range(1, self.nLtheta + 1):
                for i in range(self.nLr):
                    node_id = k * NLrtheta + j * (self.nLr + 1) + i + 1
                    eL[7].append(node_id)
        
        # Transpose to get elements (each row is one element)
        self.elemL = np.array(eL).T
        self.nemL = len(self.elemL)
    
    def generate(self):
        """
        Write local mesh files
        """
        logger.info("Writing local mesh files")
        
        from utils.format_output import format_real
        
        # 1. Write elem.l.dat
        with open('elem.l.dat', 'w') as f:
            f.write(f"{self.nemL} 8\n")
            for i, elem in enumerate(self.elemL):
                f.write(f"{i+1}")
                for node in elem:
                    f.write(f" {node}")
                f.write("\n")
        
        # 2. Write node.l.dat
        with open('node.l.dat', 'w') as f:
            f.write(f"{self.nnmL}\n")
            for i, node in enumerate(self.nodeL):
                f.write(f"{i+1}\t{format_real(node[0])}\t{format_real(node[1])}\t{format_real(node[2])}\n")
        
        logger.info("Local mesh files written successfully")
        logger.info(f"nnmG={sp.nPtsX * sp.nPtsY * sp.nPtsZ}, nnmL={self.nnmL}, ndf={3 * (sp.nPtsX * sp.nPtsY * sp.nPtsZ + self.nnmL)}")
