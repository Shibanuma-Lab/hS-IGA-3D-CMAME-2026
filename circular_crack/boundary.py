"""
Boundary conditions for circular crack problem
Direct translation from boundary[] functions in function.txt
"""

import numpy as np
import os
import sys
from const import simulation_params as sp, material_property as mp, const_local_mesh as clm
from utils.logger import logger

# Add parent directory to path for fem_data_loader import
_current_dir = os.path.dirname(os.path.abspath(__file__))
if _current_dir not in sys.path:
    sys.path.insert(0, _current_dir)

from fem_data_loader import FEMDataLoader


class Boundary:
    """
    Boundary condition handler
    Based on boundarynbc[] and boundaryebc functions from Mathematica
    """
    
    def __init__(self, local_mesh, global_mesh):
        self.local_mesh = local_mesh
        self.global_mesh = global_mesh
        
        # Boundary condition data
        self.bcG = None
        self.bcL = None
        self.load = None
        
        # Parameters
        self.nbcebc = sp.nbcebc  # 0: nodal force, 1: essential BC
        
    def define_boundary(self, local_mesh, global_mesh):
        """
        Define boundary conditions based on nbcebc setting
        """
        if self.nbcebc == 0:
            self._boundary_nbc(local_mesh.step)
        else:
            self._boundary_ebc_G(local_mesh.step)
            self._boundary_ebc_L(local_mesh.step)
        
        logger.info(f"Boundary conditions defined")
    
    def _boundary_nbc(self, step):
        """
        Boundary conditions with nodal forces (nbcebc=0)
        Translation of boundarynbc[step]
        """
        logger.info("Setting up nodal force boundary conditions")
        
        nPtsX = sp.nPtsX
        nPtsY = sp.nPtsY
        nPtsZ = sp.nPtsZ
        
        # Helper to find nodes
        def node_index_3d(i, j, k):
            """Convert 3D indices to global node number (1-based)"""
            return k * nPtsX * nPtsY + j * nPtsX + i + 1
        
        # === Global mesh boundary conditions ===
        
        # 1. bcGx0: symmetry in x-direction (x=0 plane)
        bcGx0 = []
        for k in range(nPtsZ):
            for j in range(nPtsY):
                bcGx0.append(node_index_3d(0, j, k))
        
        # 2. bcGy0: symmetry in y-direction (y=0 plane)
        bcGy0 = []
        for k in range(nPtsZ):
            for i in range(nPtsX):
                bcGy0.append(node_index_3d(i, 0, k))
        
        # 3. bcGz0: symmetry in z-direction (z=0 plane, outside crack)
        # This requires checking which elements are outside the crack
        # Simplified: all nodes at z=0
        bcGz0 = []
        for j in range(nPtsY):
            for i in range(nPtsX):
                node_id = node_index_3d(i, j, 0)
                # TODO: Filter nodes outside crack region
                bcGz0.append(node_id)
        
        # Assemble global BCs
        nbcG = len(bcGx0) + len(bcGy0) + len(bcGz0)
        
        bcG_list = [[nbcG]]
        for node in bcGx0:
            bcG_list.append([node, 1, 0.0])  # ux = 0
        for node in bcGy0:
            bcG_list.append([node, 2, 0.0])  # uy = 0
        for node in bcGz0:
            bcG_list.append([node, 3, 0.0])  # uz = 0
        
        self.bcG = bcG_list
        
        # === Applied load ===
        # Load on top surface (z = max)
        sigma_app = mp.SigmaInfinity
        
        # Simplified load distribution
        load_nodes = []
        for j in range(nPtsY):
            for i in range(nPtsX):
                node_id = node_index_3d(i, j, nPtsZ - 1)
                load_nodes.append(node_id)
        
        nloadz = len(load_nodes)
        # Distribute load equally (simplified)
        load_per_node = sigma_app * sp.WidthG * sp.HeightG / nloadz
        
        load_list = [[nloadz]]
        for node in load_nodes:
            load_list.append([node, 3, load_per_node])  # Fz
        
        self.load = load_list
        
        # === Local mesh boundary conditions ===
        
        nLr = self.local_mesh.nLr
        nLtheta = self.local_mesh.nLtheta
        HL = self.local_mesh.HL
        
        def node_index_local(i, j, k):
            """Convert local 3D indices to node number (1-based)"""
            return k * (nLr + 1) * (nLtheta + 1) + j * (nLr + 1) + i + 1
        
        # bcLx0: symmetry in y-direction (theta = 90 degrees)
        bcLx0 = []
        for k in range(HL + 1):
            for i in range(nLr + 1):
                bcLx0.append(node_index_local(i, nLtheta, k))
        
        # bcLy0: symmetry in x-direction (theta = 0)
        bcLy0 = []
        for k in range(HL + 1):
            for i in range(nLr + 1):
                bcLy0.append(node_index_local(i, 0, k))
        
        # bcLz0: symmetry in z-direction (z = 0)
        # For step <= aL, include entire bottom plane
        bcLz0 = []
        if step <= clm.aL:
            for j in range(nLtheta + 1):
                for i in range(nLr + 1):
                    bcLz0.append(node_index_local(i, j, 0))
        else:
            # No bcLz0 for propagating crack
            pass
        
        # Fixed boundaries at interfaces
        bcLr0 = []  # Inner radius (r = 0)
        for k in range(HL + 1):
            for j in range(nLtheta + 1):
                bcLr0.append(node_index_local(0, j, k))
        
        bcLr1 = []  # Outer radius (r = max)
        for k in range(HL + 1):
            for j in range(nLtheta + 1):
                bcLr1.append(node_index_local(nLr, j, k))
        
        bcLz1 = []  # Top surface (z = max)
        for j in range(nLtheta + 1):
            for i in range(nLr + 1):
                bcLz1.append(node_index_local(i, j, HL))
        
        bcLfix = list(set(bcLr0 + bcLr1 + bcLz1))
        
        # Filter symmetry BCs to remove fixed boundaries
        # Mathematica: bcLx0 = Complement[bcLx0, bcLfix];
        bcLx0 = [n for n in bcLx0 if n not in bcLfix]
        bcLy0 = [n for n in bcLy0 if n not in bcLfix]
        bcLz0 = [n for n in bcLz0 if n not in bcLfix]
        
        # Assemble local BCs
        if step == 0:
            # First step: special treatment
            # bcLr02: bottom side of inner radius (k = 0)
            bcLr02 = []
            for j in range(nLtheta + 1):
                bcLr02.append(node_index_local(0, j, 0))
            
            # bcLr01: Inner radius except top side (k < HL)
            bcLr01 = []
            for k in range(HL):
                for j in range(nLtheta + 1):
                    bcLr01.append(node_index_local(0, j, k))
            
            bcLfix0 = [n for n in bcLfix if n not in bcLr01]
            
            nbcL = (len(bcLx0) + len(bcLy0) + len(bcLz0) + len(bcLr02) +
                    len(bcLr01) * 2 + len(bcLfix0) * 3)
            
            bcL_list = [[nbcL]]
            for node in bcLx0:
                bcL_list.append([node, 1, 0.0])
            for node in bcLy0:
                bcL_list.append([node, 2, 0.0])
            for node in bcLz0:
                bcL_list.append([node, 3, 0.0])
            for node in bcLr02:
                bcL_list.append([node, 3, 0.0])
            for node in bcLr01:
                bcL_list.append([node, 1, 0.0])
                bcL_list.append([node, 2, 0.0])
            for node in bcLfix0:
                bcL_list.append([node, 1, 0.0])
                bcL_list.append([node, 2, 0.0])
                bcL_list.append([node, 3, 0.0])
        else:
            # Subsequent steps: check bcGz0[[1]] == 1
            # If crack tip is at origin (bcGz0[0] == 1), use bcLr01 separation
            # Otherwise, fix all bcLfix nodes
            if len(bcGz0) > 0 and bcGz0[0] == 1:
                # bcLr01: Inner radius except top side (k < HL)
                bcLr01 = []
                for k in range(HL):
                    for j in range(nLtheta + 1):
                        bcLr01.append(node_index_local(0, j, k))
                
                bcLfix0 = [n for n in bcLfix if n not in bcLr01]
                
                nbcL = (len(bcLx0) + len(bcLy0) + len(bcLz0) + 
                        len(bcLr01) * 2 + len(bcLfix0) * 3)
                
                bcL_list = [[nbcL]]
                for node in bcLx0:
                    bcL_list.append([node, 1, 0.0])
                for node in bcLy0:
                    bcL_list.append([node, 2, 0.0])
                for node in bcLz0:
                    bcL_list.append([node, 3, 0.0])
                for node in bcLr01:
                    bcL_list.append([node, 1, 0.0])
                    bcL_list.append([node, 2, 0.0])
                for node in bcLfix0:
                    bcL_list.append([node, 1, 0.0])
                    bcL_list.append([node, 2, 0.0])
                    bcL_list.append([node, 3, 0.0])
            else:
                # Fix all bcLfix nodes in all directions
                nbcL = len(bcLx0) + len(bcLy0) + len(bcLz0) + len(bcLfix) * 3
                
                bcL_list = [[nbcL]]
                for node in bcLx0:
                    bcL_list.append([node, 1, 0.0])
                for node in bcLy0:
                    bcL_list.append([node, 2, 0.0])
                for node in bcLz0:
                    bcL_list.append([node, 3, 0.0])
                for node in bcLfix:
                    bcL_list.append([node, 1, 0.0])
                    bcL_list.append([node, 2, 0.0])
                    bcL_list.append([node, 3, 0.0])
        
        self.bcL = bcL_list
    
    def _boundary_ebc_G(self, step):
        """
        Essential boundary conditions for global mesh (nbcebc=1)
        Translation of boundaryebcG[step]
        
        Uses visualization mesh (node_visual, elem_visual) to determine crack region,
        then applies BCs to control points (nodeG)
        """
        logger.info(f"Setting up essential boundary conditions for global mesh (step={step})")
        
        nodeG = self.global_mesh.nodeG
        node_visual = self.global_mesh.node_visual
        elem_visual = self.global_mesh.elem_visual
        elemG = self.global_mesh.elemG
        
        nPtsX = sp.nPtsX
        nPtsY = sp.nPtsY
        nPtsZ = sp.nPtsZ
        
        # Find nodes on boundaries (control points)
        tol = 1e-10
        
        # bcGx0: x = 0
        bcGx0 = np.where(np.abs(nodeG[:, 0]) < tol)[0] + 1  # 1-based
        
        # bcGy0: y = 0
        bcGy0 = np.where(np.abs(nodeG[:, 1]) < tol)[0] + 1
        
        # bcGz0: z = 0 outside crack region
        # Use visualization mesh to find elements outside crack
        # eGout: elements where max radius of visual nodes > step*hL
        nelemU = len(np.unique(self.global_mesh.uKnot)) - 1
        nelemV = len(np.unique(self.global_mesh.vKnot)) - 1
        
        eGout = []
        for e in range(nelemU * nelemV):
            if e >= len(elem_visual):
                break
            # Get visual nodes of this element
            elem_nodes = elem_visual[e]
            # Calculate max radius in xy-plane
            max_r = 0.0
            for node_id in elem_nodes:
                if node_id > 0 and node_id <= len(node_visual):
                    node_pos = node_visual[node_id - 1]
                    r = np.sqrt(node_pos[0]**2 + node_pos[1]**2)
                    max_r = max(max_r, r)
            
            if max_r > step * clm.hL:
                eGout.append(e)
        
        # nGout: control points in elements outside crack
        nGout = set()
        for e in eGout:
            if e < len(elemG):
                nGout.update(elemG[e])
        nGout = sorted(list(nGout))
        
        # bcGz0: intersection of nGout with z=0 plane control points
        z0_nodes = set(range(1, nPtsX * nPtsY + 1))
        bcGz0 = sorted(list(z0_nodes.intersection(set(nGout))))
        
        # Outer boundaries (for prescribed displacement)
        # bcGx1: x = max (Mathematica: (nPtsX)*Range[(nPtsX)*(nPtsZ)])
        bcGx1 = [(m + 1) * nPtsX for m in range(nPtsX * nPtsZ)]
        
        # bcGy1: y = max (Mathematica: Flatten[Map[# - Range[0, nPtsX] &, (nPtsX)^2*Range[nPtsZ]]])
        bcGy1 = []
        for k in range(1, nPtsZ + 1):
            base = k * nPtsX * nPtsX
            for i in range(nPtsX + 1):
                bcGy1.append(base - i)
        
        # bcGz1: z = max (Mathematica: (nPtsX)^2*(nPtsZ) - Range[(nPtsX)^2] + 1)
        # Note: Mathematica Range[n] = [1, 2, ..., n], so we use range(1, n+1)
        bcGz1 = [nPtsX * nPtsX * nPtsZ - m + 1 for m in range(1, nPtsX * nPtsX + 1)]
        
        bcGs = sorted(list(set(bcGx1 + bcGy1 + bcGz1)))
        
        # Remove outer boundary nodes from inner symmetry BC
        bcGx0 = [n for n in bcGx0 if n not in bcGs]
        bcGy0 = [n for n in bcGy0 if n not in bcGs]
        bcGz0 = [n for n in bcGz0 if n not in bcGs]
        
        # Get displacement values for outer boundary nodes
        bc_values = self._get_fem_bc_values(bcGs, step)
        
        # Assemble BCs
        nbcG = len(bcGx0) + len(bcGy0) + len(bcGz0) + 3 * len(bcGs)
        
        bcG_list = [[nbcG]]
        for node in bcGx0:
            bcG_list.append([int(node), 1, 0.0])
        for node in bcGy0:
            bcG_list.append([int(node), 2, 0.0])
        for node in bcGz0:
            bcG_list.append([int(node), 3, 0.0])
        
        # For nodes on outer boundary, use FEM-interpolated values
        for node in bcGs:
            node_idx = int(node) - 1  # Convert to 0-based
            if node_idx < len(bc_values):
                bcG_list.append([int(node), 1, bc_values[node_idx, 0]])
                bcG_list.append([int(node), 2, bc_values[node_idx, 1]])
                bcG_list.append([int(node), 3, bc_values[node_idx, 2]])
            else:
                bcG_list.append([int(node), 1, 0.0])
                bcG_list.append([int(node), 2, 0.0])
                bcG_list.append([int(node), 3, 0.0])
        
        self.bcG = bcG_list
        
        # No applied load for essential BC case
        self.load = [[0]]
    
    def _get_fem_bc_values(self, node_indices, step):
        """
        Get boundary condition values from FEM interpolation
        
        Args:
            node_indices: Array of node indices (1-based)
            step: Time step number
        
        Returns:
            Array of shape (n_total_nodes, 3) with displacement values
            
        Note: In Mathematica code, FEM data is ALWAYS used for bcGs nodes,
        regardless of step value. The step parameter only determines which
        time slice of FEM data to use.
        """
        n_nodes = len(self.global_mesh.nodeG)
        bc_values = np.zeros((n_nodes, 3))
        
        # Try to load FEM data
        fem_loader = FEMDataLoader()
        
        if not fem_loader.load_data(sp.V):
            logger.warning(f"Could not load FEM data for V={sp.V}, using zero displacement")
            return bc_values
        
        # Get displacement at this step
        # Note: step is used as index into FEM data (step+1 in Mathematica, but we use 0-based)
        disp_2d = fem_loader.get_displacement_at_step(step)
        if disp_2d is None:
            logger.warning(f"No FEM displacement data for step {step}, using zero")
            return bc_values
        
        # Interpolate 2D FEM values to 3D IGA nodes
        disp_3d = fem_loader.interpolate_2d_to_3d(disp_2d, self.global_mesh.nodeG)
        
        return disp_3d
    
    def _boundary_ebc_L(self, step):
        """
        Essential boundary conditions for local mesh (nbcebc=1)
        Translation of boundaryebcL[step]
        """
        logger.info("Setting up essential boundary conditions for local mesh")
        
        nLr = self.local_mesh.nLr
        nLtheta = self.local_mesh.nLtheta
        HL = self.local_mesh.HL
        
        def node_index_local(i, j, k):
            return k * (nLr + 1) * (nLtheta + 1) + j * (nLr + 1) + i + 1
        
        # Same as in _boundary_nbc
        bcLx0 = []
        for k in range(HL + 1):
            for i in range(nLr + 1):
                bcLx0.append(node_index_local(i, nLtheta, k))
        
        bcLy0 = []
        for k in range(HL + 1):
            for i in range(nLr + 1):
                bcLy0.append(node_index_local(i, 0, k))
        
        # bcLz0: symmetry in z-direction (z = 0)
        # Mathematica: Range[step + 1, nLr + 1] means i starts from step (in 0-based)
        bcLz0 = []
        if step <= clm.aL:
            i_start = step  # For step=0: i>=0, step=2: i>=2
            for j in range(nLtheta + 1):
                for i in range(i_start, nLr + 1):
                    bcLz0.append(node_index_local(i, j, 0))
        else:
            i_start = clm.aL
            for j in range(nLtheta + 1):
                for i in range(i_start, nLr + 1):
                    bcLz0.append(node_index_local(i, j, 0))
        
        bcLr0 = []
        for k in range(HL + 1):
            for j in range(nLtheta + 1):
                bcLr0.append(node_index_local(0, j, k))
        
        bcLr1 = []
        for k in range(HL + 1):
            for j in range(nLtheta + 1):
                bcLr1.append(node_index_local(nLr, j, k))
        
        bcLz1 = []
        for j in range(nLtheta + 1):
            for i in range(nLr + 1):
                bcLz1.append(node_index_local(i, j, HL))
        
        bcLfix = list(set(bcLr0 + bcLr1 + bcLz1))
        
        # Filter symmetry BCs (Mathematica: bcLx0 = Complement[bcLx0, bcLfix])
        bcLx0 = [n for n in bcLx0 if n not in bcLfix]
        bcLy0 = [n for n in bcLy0 if n not in bcLfix]
        bcLz0 = [n for n in bcLz0 if n not in bcLfix]
        
        if step == 0:
            # bcLr02: bottom side of inner radius (k = 0)
            bcLr02 = []
            for j in range(nLtheta + 1):
                bcLr02.append(node_index_local(0, j, 0))
            
            # bcLr01: inner side except for top side (k < HL)
            bcLr01 = []
            for k in range(HL):
                for j in range(nLtheta + 1):
                    bcLr01.append(node_index_local(0, j, k))
            
            bcLfix0 = [n for n in bcLfix if n not in bcLr01]
            
            nbcL = (len(bcLx0) + len(bcLy0) + len(bcLz0) + len(bcLr02) +
                    len(bcLr01) * 2 + len(bcLfix0) * 3)
            
            bcL_list = [[nbcL]]
            for node in bcLx0:
                bcL_list.append([node, 1, 0.0])
            for node in bcLy0:
                bcL_list.append([node, 2, 0.0])
            for node in bcLz0:
                bcL_list.append([node, 3, 0.0])
            for node in bcLr02:
                bcL_list.append([node, 3, 0.0])
            for node in bcLr01:
                bcL_list.append([node, 1, 0.0])
                bcL_list.append([node, 2, 0.0])
            for node in bcLfix0:
                bcL_list.append([node, 1, 0.0])
                bcL_list.append([node, 2, 0.0])
                bcL_list.append([node, 3, 0.0])
        else:
            # Subsequent steps: check if bcGz0[[1]] == 1 OR step <= aL
            # This matches Mathematica: If[bcGz0[[1]] == 1 || step <= aL, ...]
            # Get bcGz0 from global mesh BC calculation
            nodeG = self.global_mesh.nodeG
            tol = 1e-10
            all_z0 = np.where(np.abs(nodeG[:, 2]) < tol)[0]
            crack_radius = step * clm.hL
            bcGz0 = []
            for idx in all_z0:
                r = np.sqrt(nodeG[idx, 0]**2 + nodeG[idx, 1]**2)
                if r > crack_radius:
                    bcGz0.append(idx + 1)
            
            if (len(bcGz0) > 0 and bcGz0[0] == 1) or (step <= clm.aL):
                # bcLr01: Inner radius except top side (k < HL)
                bcLr01 = []
                for k in range(HL):
                    for j in range(nLtheta + 1):
                        bcLr01.append(node_index_local(0, j, k))
                
                bcLfix0 = [n for n in bcLfix if n not in bcLr01]
                
                nbcL = (len(bcLx0) + len(bcLy0) + len(bcLz0) + 
                        len(bcLr01) * 2 + len(bcLfix0) * 3)
                
                bcL_list = [[nbcL]]
                for node in bcLx0:
                    bcL_list.append([node, 1, 0.0])
                for node in bcLy0:
                    bcL_list.append([node, 2, 0.0])
                for node in bcLz0:
                    bcL_list.append([node, 3, 0.0])
                for node in bcLr01:
                    bcL_list.append([node, 1, 0.0])
                    bcL_list.append([node, 2, 0.0])
                for node in bcLfix0:
                    bcL_list.append([node, 1, 0.0])
                    bcL_list.append([node, 2, 0.0])
                    bcL_list.append([node, 3, 0.0])
            else:
                # Fix all bcLfix nodes in all directions
                nbcL = len(bcLx0) + len(bcLy0) + len(bcLz0) + len(bcLfix) * 3
                
                bcL_list = [[nbcL]]
                for node in bcLx0:
                    bcL_list.append([node, 1, 0.0])
                for node in bcLy0:
                    bcL_list.append([node, 2, 0.0])
                for node in bcLz0:
                    bcL_list.append([node, 3, 0.0])
                for node in bcLfix:
                    bcL_list.append([node, 1, 0.0])
                    bcL_list.append([node, 2, 0.0])
                    bcL_list.append([node, 3, 0.0])
        
        self.bcL = bcL_list
    
    def generate(self):
        """
        Write boundary condition files
        """
        logger.info("Writing boundary condition files")
        
        from utils.format_output import format_bc_line, format_real
        
        # Write bc.g.dat
        with open('bc.g.dat', 'w') as f:
            for bc in self.bcG:
                if len(bc) == 1:
                    f.write(f"{bc[0]}\n")
                else:
                    f.write(format_bc_line(bc[0], bc[1], bc[2]) + "\n")
        
        # Write bc.l.dat
        with open('bc.l.dat', 'w') as f:
            for bc in self.bcL:
                if len(bc) == 1:
                    f.write(f"{bc[0]}\n")
                else:
                    f.write(format_bc_line(bc[0], bc[1], bc[2]) + "\n")
        
        # Write load.dat
        with open('load.dat', 'w') as f:
            for load in self.load:
                if len(load) == 1:
                    f.write(f"{load[0]}\n")
                else:
                    f.write(f"{load[0]}\t{load[1]}\t{format_real(load[2])}\n")
        
        logger.info("Boundary condition files written successfully")
