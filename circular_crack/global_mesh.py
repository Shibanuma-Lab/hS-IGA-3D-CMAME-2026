"""
Global mesh generation using IGA (Isogeometric Analysis)
Direct translation from CircularCrackIGA[] in function.txt
"""

import numpy as np
from const import const_global_mesh as cgm, simulation_params as sp
from utils.logger import logger
from utils.nurbs_utils import find_span, basis_funs


class GlobalMesh:
    """
    Global mesh generator using B-spline/NURBS basis
    Based on CircularCrackIGA[] from Mathematica
    """
    
    def __init__(self, step):
        self.step = step
        
        # B-spline degrees
        self.p = cgm.p
        self.q = cgm.q
        self.r = cgm.r
        
        # Number of control points
        self.nPtsX = sp.nPtsX
        self.nPtsY = sp.nPtsY
        self.nPtsZ = sp.nPtsZ
        
        # Element size
        if sp.static_mode:
            # Static mode: use non-uniform spacing
            self.hGx = sp.hGX
            self.hGy = sp.hGY
            self.hGz = sp.hGZ
            self.hG = sp.hG  # Keep for compatibility
        else:
            # Dynamic mode: uniform spacing
            self.hG = sp.hG
        
        # Mesh data
        self.nodeG = None
        self.elemG = None  # element connectivity
        self.weights = None
        self.indexG = None
        
        # Knot vectors
        self.uKnot = None
        self.vKnot = None
        self.wKnot = None
        
        # For visualization
        self.node_visual = None
        self.elem_visual = None
        
    def make_global_mesh(self):
        """
        Generate global IGA mesh
        Translation of CircularCrackIGA[hG]
        """
        logger.info("Generating global IGA mesh")
        
        # Helper function: set very small values to zero
        def ad(x):
            return 0.0 if abs(x) < 1e-9 * sp.WidthG else x
        
        # 1. Generate nodal coordinates (control points)
        if sp.static_mode:
            # Static mode: non-uniform element spacing
            nodeGx = self.hGx * np.arange(self.nPtsX)
            nodeGy = self.hGy * np.arange(self.nPtsY)
            nodeGz = self.hGz * np.arange(self.nPtsZ)
        else:
            # Dynamic mode: uniform element spacing
            nodeGx = self.hG * np.arange(self.nPtsX)
            nodeGy = self.hG * np.arange(self.nPtsY)
            nodeGz = self.hG * np.arange(self.nPtsZ)
        
        # Create 2D grid for XY plane
        nodeGxy = []
        for y in nodeGy:
            for x in nodeGx:
                nodeGxy.append([x, y])
        nodeGxy = np.array(nodeGxy)
        
        # Create 3D control points
        nodeG_list = []
        for z in nodeGz:
            for xy in nodeGxy:
                nodeG_list.append([ad(xy[0]), ad(xy[1]), ad(z)])
        
        self.nodeG = np.array(nodeG_list)
        nnodeG = len(self.nodeG)
        
        logger.info(f"Control points generated: {nnodeG} nodes")
        
        # 2. Generate knot vectors
        # For open B-spline: repeat first and last knots (p+1) times
        knotUTemp = np.linspace(0, 1, self.nPtsX - self.p + 1)
        knotVTemp = np.linspace(0, 1, self.nPtsY - self.q + 1)
        knotWTemp = np.linspace(0, 1, self.nPtsZ - self.r + 1)
        
        # Add repeated knots at boundaries
        self.uKnot = np.concatenate([[0, 0], knotUTemp, [1, 1]])
        self.vKnot = np.concatenate([[0, 0], knotVTemp, [1, 1]])
        self.wKnot = np.concatenate([[0, 0], knotWTemp, [1, 1]])
        
        logger.info(f"Knot vectors: U={len(self.uKnot)}, V={len(self.vKnot)}, W={len(self.wKnot)}")
        
        # 3. Generate weights (all 1.0 for B-spline)
        self.weights = np.ones(self.nPtsX * self.nPtsY * self.nPtsZ)
        
        # 4. Generate element connectivity
        self._generate_elements()
        
        # 5. Generate visualization mesh
        self._build_visual_mesh()
        
        logger.info(f"Global mesh complete: {len(self.elemG)} elements")
    
    def _generate_elements(self):
        """
        Generate element connectivity using FindSpan
        Translation from Mathematica code section
        """
        noU = self.nPtsX
        noV = self.nPtsY
        noW = self.nPtsZ
        
        p, q, r = self.p, self.q, self.r
        
        # 1. Unique knots & number of elements
        uniqU = np.unique(self.uKnot)
        uniqV = np.unique(self.vKnot)
        uniqW = np.unique(self.wKnot)
        
        nelemU = len(uniqU) - 1
        nelemV = len(uniqV) - 1
        nelemW = len(uniqW) - 1
        
        # 2. Parameter ranges for each element
        elRangeU = [[uniqU[i], uniqU[i+1]] for i in range(nelemU)]
        elRangeV = [[uniqV[j], uniqV[j+1]] for j in range(nelemV)]
        elRangeW = [[uniqW[k], uniqW[k+1]] for k in range(nelemW)]
        
        # 3. Connectivity via FindSpan
        elConnU = []
        for i in range(nelemU):
            xi_mid = 0.5 * (elRangeU[i][0] + elRangeU[i][1])
            span = find_span(noU - 1, p, xi_mid, self.uKnot)
            elConnU.append(list(range(span - p + 1, span + 2)))
        
        elConnV = []
        for j in range(nelemV):
            eta_mid = 0.5 * (elRangeV[j][0] + elRangeV[j][1])
            span = find_span(noV - 1, q, eta_mid, self.vKnot)
            elConnV.append(list(range(span - q + 1, span + 2)))
        
        elConnW = []
        for k in range(nelemW):
            zeta_mid = 0.5 * (elRangeW[k][0] + elRangeW[k][1])
            span = find_span(noW - 1, r, zeta_mid, self.wKnot)
            elConnW.append(list(range(span - r + 1, span + 2)))
        
        # 4. Global node numbering chan[w,v,u]
        # Mathematica: chan = ArrayReshape[Range[noU*noV*noW], {noW, noV, noU}]
        chan = np.arange(1, noU * noV * noW + 1).reshape(noW, noV, noU)
        
        # 5. Build element list
        elements = []
        for w in range(nelemW):
            for v in range(nelemV):
                for u in range(nelemU):
                    elem_nodes = []
                    for wk in elConnW[w]:
                        for vj in elConnV[v]:
                            for ui in elConnU[u]:
                                # Convert to 0-based index for chan access
                                node_id = chan[wk-1, vj-1, ui-1]
                                elem_nodes.append(node_id)
                    elements.append(elem_nodes)
        
        self.elemG = np.array(elements, dtype=int)
        self.nelem = len(self.elemG)
        
        # 6. Generate index (element parameter range indices)
        self.indexG = []
        for w in range(nelemW):
            for v in range(nelemV):
                for u in range(nelemU):
                    self.indexG.append([u+1, v+1, w+1])  # 1-based for Fortran
        self.indexG = np.array(self.indexG, dtype=int)
    
    def _build_visual_mesh(self):
        """
        Build visualization mesh for VTU output
        Translation of buildVisual3D[] from Mathematica
        Uses NURBS surface evaluation at knot locations
        """
        # Get unique knot values
        uKnotVec = np.unique(self.uKnot)
        vKnotVec = np.unique(self.vKnot)
        wKnotVec = np.unique(self.wKnot)
        
        noKnotsU = len(uKnotVec)
        noKnotsV = len(vKnotVec)
        noKnotsW = len(wKnotVec)
        
        # Generate visualization nodes using NURBS surface evaluation
        node_visual = []
        for wk in range(noKnotsW):
            for vj in range(noKnotsV):
                for ui in range(noKnotsU):
                    u = uKnotVec[ui]
                    v = vKnotVec[vj]
                    w = wKnotVec[wk]
                    
                    # Evaluate NURBS solid at parameter point (u, v, w)
                    point = self._evaluate_nurbs_point(u, v, w)
                    node_visual.append(point)
        
        self.node_visual = np.array(node_visual)
        
        # Generate hex elements for visualization
        nx = noKnotsU - 1
        ny = noKnotsV - 1
        nz = noKnotsW - 1
        nnodexy = noKnotsU * noKnotsV
        
        elem_visual = []
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    # 8-node hexahedron connectivity
                    n1 = k * nnodexy + j * noKnotsU + i + 1
                    n2 = n1 + 1
                    n3 = n1 + noKnotsU + 1
                    n4 = n1 + noKnotsU
                    n5 = n1 + nnodexy
                    n6 = n5 + 1
                    n7 = n5 + noKnotsU + 1
                    n8 = n5 + noKnotsU
                    
                    elem_visual.append([n1, n2, n3, n4, n5, n6, n7, n8])
        
        self.elem_visual = np.array(elem_visual, dtype=int)
    
    def _evaluate_nurbs_point(self, u, v, w):
        """
        Evaluate NURBS solid at parameter coordinates (u, v, w)
        Translation of SolidPoint[] from Mathematica
        
        Returns: [x, y, z] physical coordinates
        """
        # Find spans
        uspan = find_span(self.nPtsX - 1, self.p, u, self.uKnot)
        vspan = find_span(self.nPtsY - 1, self.q, v, self.vKnot)
        wspan = find_span(self.nPtsZ - 1, self.r, w, self.wKnot)
        
        # Compute basis functions
        Nu = basis_funs(uspan, u, self.p, self.uKnot)
        Nv = basis_funs(vspan, v, self.q, self.vKnot)
        Nw = basis_funs(wspan, w, self.r, self.wKnot)
        
        # Local control point indices
        uind = uspan - self.p
        vind = vspan - self.q
        wind = wspan - self.r
        
        # Compute weighted sum
        point = np.zeros(3)
        weight_sum = 0.0
        
        for k in range(self.r + 1):
            for j in range(self.q + 1):
                for i in range(self.p + 1):
                    # Global control point index (0-based)
                    idx = (wind + k) * self.nPtsX * self.nPtsY + (vind + j) * self.nPtsX + (uind + i)
                    
                    # Basis function value
                    N = Nu[i] * Nv[j] * Nw[k]
                    
                    # Weight (all 1.0 for B-spline)
                    w_cp = self.weights[idx]
                    
                    # Add weighted contribution
                    point += N * w_cp * self.nodeG[idx]
                    weight_sum += N * w_cp
        
        # Normalize by weight sum
        point /= weight_sum
        
        return point
    
    def generate(self):
        """
        Write mesh files in the format expected by SFEM solver
        """
        logger.info("Writing global mesh files")
        
        from utils.format_output import format_real
        
        # 1. Write elem.g.dat
        nelem = len(self.elemG)
        nnodes_per_elem = self.elemG.shape[1]
        
        with open('elem.g.dat', 'w') as f:
            f.write(f"{nelem} {nnodes_per_elem}\n")
            for i, elem in enumerate(self.elemG):
                f.write(f"{i+1}")
                for node in elem:
                    f.write(f" {node}")
                f.write("\n")
        
        # 2. Write index.g.dat
        with open('index.g.dat', 'w') as f:
            f.write(f"{len(self.indexG)}\n")
            for idx in self.indexG:
                f.write(f"{idx[0]} {idx[1]} {idx[2]}\n")
        
        # 3. Write weights.g.dat (no node IDs, only values)
        with open('weights.g.dat', 'w') as f:
            f.write(f"{len(self.weights)}\n")
            for w in self.weights:
                f.write(f"{format_real(w)}\n")
        
        # 4. Write node.g.dat
        with open('node.g.dat', 'w') as f:
            f.write(f"{len(self.nodeG)}\n")
            for i, node in enumerate(self.nodeG):
                f.write(f"{i+1}\t{format_real(node[0])}\t{format_real(node[1])}\t{format_real(node[2])}\n")
        
        # 5. Write visualization files
        with open('elem.v.dat', 'w') as f:
            f.write(f"{len(self.elem_visual)} 8\n")
            for i, elem in enumerate(self.elem_visual):
                f.write(f"{i+1}")
                for node in elem:
                    f.write(f" {node}")
                f.write("\n")
        
        with open('node.v.dat', 'w') as f:
            f.write(f"{len(self.node_visual)}\n")
            for i, node in enumerate(self.node_visual):
                f.write(f"{i+1}\t{format_real(node[0])}\t{format_real(node[1])}\t{format_real(node[2])}\n")
        
        logger.info("Global mesh files written successfully")
