"""
Complete J-integral calculation matching Mathematica implementation
Includes mesh generation and data extension from 0-90° to 0-180°
"""

import numpy as np
import csv
from pathlib import Path
from numpy.polynomial.legendre import leggauss
from const import material_property as mp
from const import const_local_mesh as clm
from const import simulation_params as sp

Path("logs").mkdir(exist_ok=True)

try:
    from fem_data_loader import FEMDataLoader
except ImportError:  # pragma: no cover - package import fallback
    from .fem_data_loader import FEMDataLoader


class JIntegralCalculator:
    """Direct translation of Mathematica Jintegral[] function with mesh generation"""
    
    def __init__(
        self,
        step_start=1,
        step_end=1,
        Rj0=1.5,
        Rj1=1.515,
        Wj0=1.0,
        Wj1=1.01,
        v=1000.0,
        result_root=None,
        sigma_app=None,
        c=None,
        nu=None,
        EE=None,
        rho=None,
        aL=None,
        lL=None,
        d_theta=None,
        hL=None,
        HL=None,
        ngp=None,
    ):
        """
        Initialize J-integral calculator with configurable parameters.
        
        Args:
            step_start: Starting step number
            step_end: Ending step number
            Rj0: Inner radius for J-integral domain (normalized by hL)
            Rj1: Outer radius for J-integral domain
            Wj0: Inner width parameter
            Wj1: Outer width parameter
            v: Crack velocity (m/s)
        """
        # Import parameters from const modules
        self.sigma_app = mp.SigmaInfinity if sigma_app is None else sigma_app
        self.c = sp.c if c is None else c
        self.nu = mp.Nu if nu is None else nu
        self.EE = mp.EE if EE is None else EE
        self.rho = mp.Rho if rho is None else rho

        self.aL = clm.aL if aL is None else aL
        self.lL = clm.lL if lL is None else lL
        self.d_theta = clm.d_theta if d_theta is None else d_theta
        self.hL = clm.hL if hL is None else hL
        self.HL = clm.HL if HL is None else HL  # Height parameter (厚度方向单元数)
        self.result_root = Path("results") if result_root is None else Path(result_root)
        
        # Configurable J-integral parameters
        self.Rj0 = Rj0
        self.Rj1 = Rj1
        self.Wj0 = Wj0
        self.Wj1 = Wj1
        
        self.v = v
        self.ngp = 2 if ngp is None else int(ngp)
        
        # Step range
        self.stepini = step_start
        self.stepend = step_end
        
        # Calculated parameters
        self.nLr = self.aL + self.lL
        self.nL_theta = round(90 / self.d_theta)
        
        # Hooke tensor
        factor = self.EE / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        self.de = factor * np.array([
            [1.0 - self.nu, self.nu, self.nu, 0, 0, 0],
            [self.nu, 1.0 - self.nu, self.nu, 0, 0, 0],
            [self.nu, self.nu, 1.0 - self.nu, 0, 0, 0],
            [0, 0, 0, 0.5 - self.nu, 0, 0],
            [0, 0, 0, 0, 0.5 - self.nu, 0],
            [0, 0, 0, 0, 0, 0.5 - self.nu]
        ])
        
        # Gauss points and weights
        self.setup_gauss()
        
        # Results storage
        self.JAll = []
        self.KIDAll = []
        
    def setup_gauss(self):
        """Setup Gauss integration points and weights"""
        gp, gw = leggauss(self.ngp)
        gp = np.asarray(gp, dtype=np.float64)
        gw = np.asarray(gw, dtype=np.float64)
        
        # 3D Gauss points: [psi, eta, zeta]
        # Match Mathematica order: eta, zeta, psi
        self.psi_eta_zeta = []
        self.weight = []
        for i_eta, eta in enumerate(gp):
            for i_zeta, zeta in enumerate(gp):
                for i_psi, psi in enumerate(gp):
                    self.psi_eta_zeta.append([psi, eta, zeta])
                    self.weight.append(gw[i_psi] * gw[i_eta] * gw[i_zeta])
        
        # Shape functions at Gauss points
        self.nn = [self.shp(p) for p in self.psi_eta_zeta]
        self.Dnn = [self.Dshp(p) for p in self.psi_eta_zeta]
    
    def shp(self, coords):
        """Shape functions for C3D8 element"""
        psi, eta, zeta = coords
        return (1.0/8.0) * np.array([
            (1 - psi) * (1 - eta) * (1 - zeta),
            (1 + psi) * (1 - eta) * (1 - zeta),
            (1 + psi) * (1 + eta) * (1 - zeta),
            (1 - psi) * (1 + eta) * (1 - zeta),
            (1 - psi) * (1 - eta) * (1 + zeta),
            (1 + psi) * (1 - eta) * (1 + zeta),
            (1 + psi) * (1 + eta) * (1 + zeta),
            (1 - psi) * (1 + eta) * (1 + zeta)
        ])
    
    def Dshp(self, coords):
        """Shape function derivatives for C3D8 element"""
        psi, eta, zeta = coords
        return (1.0/8.0) * np.array([
            # d/d_psi
            [-(1 - eta) * (1 - zeta), (1 - eta) * (1 - zeta), (1 + eta) * (1 - zeta), -(1 + eta) * (1 - zeta),
             -(1 - eta) * (1 + zeta), (1 - eta) * (1 + zeta), (1 + eta) * (1 + zeta), -(1 + eta) * (1 + zeta)],
            # d/d_eta
            [-(1 - psi) * (1 - zeta), -(1 + psi) * (1 - zeta), (1 + psi) * (1 - zeta), (1 - psi) * (1 - zeta),
             -(1 - psi) * (1 + zeta), -(1 + psi) * (1 + zeta), (1 + psi) * (1 + zeta), (1 - psi) * (1 + zeta)],
            # d/d_zeta
            [-(1 - psi) * (1 - eta), -(1 + psi) * (1 - eta), -(1 + psi) * (1 + eta), -(1 - psi) * (1 + eta),
             (1 - psi) * (1 - eta), (1 + psi) * (1 - eta), (1 + psi) * (1 + eta), (1 - psi) * (1 + eta)]
        ])
    
    def generate_mesh(self, step):
        """
        Generate nodeLX mesh matching Mathematica
        Mathematica: nodeLX = Flatten[Outer[{ad[#3*Cos[#2]], ad[#3*Sin[#2]], ad[#1]}&, 
                                             nodeLz, nodeLθX, nodeLr], 2]
        """
        # Generate coordinate arrays
        if step <= self.aL:
            nodeLr = self.hL * np.arange(0, self.nLr + 1)
        else:
            nodeLr = self.hL * (step + np.arange(-self.aL, self.lL + 1))
        
        nodeLz = self.hL * np.arange(0, self.HL + 1)
        nodeLtheta_X = np.arange(0, np.pi + 1e-10, (0.5 * np.pi) / self.nL_theta)
        
        # Generate mesh using outer product
        # Order: z, theta, r (matching Mathematica's Outer order)
        nodeLX = []
        for z in nodeLz:
            for theta in nodeLtheta_X:
                for r in nodeLr:
                    x = r * np.cos(theta)
                    y = r * np.sin(theta)
                    # Adjust to zero for very small values
                    if abs(x) < 1e-9:
                        x = 0.0
                    if abs(y) < 1e-9:
                        y = 0.0
                    if abs(z) < 1e-9:
                        z = 0.0
                    nodeLX.append([x, y, z])
        
        return np.array(nodeLX)

    def generate_half_mesh(self, step):
        """
        Generate the original 0-90 degree local mesh nodes.
        """
        if step <= self.aL:
            nodeLr = self.hL * np.arange(0, self.nLr + 1)
        else:
            nodeLr = self.hL * (step + np.arange(-self.aL, self.lL + 1))

        nodeLz = self.hL * np.arange(0, self.HL + 1)
        nodeLtheta = np.arange(0, 0.5 * np.pi + 1e-10, (0.5 * np.pi) / self.nL_theta)

        nodeL = []
        for z in nodeLz:
            for theta in nodeLtheta:
                for r in nodeLr:
                    x = r * np.cos(theta)
                    y = r * np.sin(theta)
                    if abs(x) < 1e-9:
                        x = 0.0
                    if abs(y) < 1e-9:
                        y = 0.0
                    if abs(z) < 1e-9:
                        z = 0.0
                    nodeL.append([x, y, z])

        return np.array(nodeL)
    
    def generate_elements(self, step):
        """
        Generate elemLX matching Mathematica EXACTLY
        Mathematica formula:
        eLX[[1]] = Flatten[Outer[#1 + #2 + #3 &, NLXrθ*Range[0, HL-1], 
                                  (nLr+1)*Range[0, 2*nLθ-1], Range[1, nLr]]]
        eLX[[2]] = ... Range[2, nLr+1]
        eLX[[3]] = ... Range[1, 2*nLθ], Range[2, nLr+1]
        eLX[[4]] = ... Range[1, 2*nLθ], Range[1, nLr]
        eLX[[5-8]] same but Range[1, HL] for z
        """
        if step <= self.aL:
            nLr = self.nLr
        else:
            nLr = self.aL + self.lL
        
        num_r = nLr + 1
        num_theta = 2 * self.nL_theta + 1
        NLXrtheta = num_r * num_theta
        
        # Generate 8 node lists (matching Mathematica's eLX[[1]] to eLX[[8]])
        eLX = [[] for _ in range(8)]
        
        # eLX[[1]]: bottom-left-near
        # Range[0, HL-1] * NLXrθ + Range[0, 2*nLθ-1] * (nLr+1) + Range[1, nLr]
        for iz in range(self.HL):
            for itheta in range(2 * self.nL_theta):
                for ir in range(1, nLr + 1):  # Range[1, nLr] in Mathematica
                    node = iz * NLXrtheta + itheta * num_r + ir
                    eLX[0].append(node)
        
        # eLX[[2]]: bottom-right-near (r+1)
        for iz in range(self.HL):
            for itheta in range(2 * self.nL_theta):
                for ir in range(2, nLr + 2):  # Range[2, nLr+1]
                    node = iz * NLXrtheta + itheta * num_r + ir
                    eLX[1].append(node)
        
        # eLX[[3]]: bottom-right-far (r+1, theta+1)
        for iz in range(self.HL):
            for itheta in range(1, 2 * self.nL_theta + 1):  # Range[1, 2*nLθ]
                for ir in range(2, nLr + 2):
                    node = iz * NLXrtheta + itheta * num_r + ir
                    eLX[2].append(node)
        
        # eLX[[4]]: bottom-left-far (theta+1)
        for iz in range(self.HL):
            for itheta in range(1, 2 * self.nL_theta + 1):
                for ir in range(1, nLr + 1):
                    node = iz * NLXrtheta + itheta * num_r + ir
                    eLX[3].append(node)
        
        # eLX[[5]]: top-left-near (z+1)
        for iz in range(1, self.HL + 1):  # Range[1, HL]
            for itheta in range(2 * self.nL_theta):
                for ir in range(1, nLr + 1):
                    node = iz * NLXrtheta + itheta * num_r + ir
                    eLX[4].append(node)
        
        # eLX[[6]]: top-right-near (z+1, r+1)
        for iz in range(1, self.HL + 1):
            for itheta in range(2 * self.nL_theta):
                for ir in range(2, nLr + 2):
                    node = iz * NLXrtheta + itheta * num_r + ir
                    eLX[5].append(node)
        
        # eLX[[7]]: top-right-far (z+1, r+1, theta+1)
        for iz in range(1, self.HL + 1):
            for itheta in range(1, 2 * self.nL_theta + 1):
                for ir in range(2, nLr + 2):
                    node = iz * NLXrtheta + itheta * num_r + ir
                    eLX[6].append(node)
        
        # eLX[[8]]: top-left-far (z+1, theta+1)
        for iz in range(1, self.HL + 1):
            for itheta in range(1, 2 * self.nL_theta + 1):
                for ir in range(1, nLr + 1):
                    node = iz * NLXrtheta + itheta * num_r + ir
                    eLX[7].append(node)
        
        # Transpose to get elements (each row is one element with 8 nodes)
        elemLX = np.array(eLX).T
        
        return elemLX
    
    def extend_data_to_full_circle(self, data_half, step):
        """
        Extend data from 0-90° to 0-180° using mirror symmetry
        Mathematica: disLGX = Flatten[MapThread[Join[#1, #2]&, {disLGX1, disLGX2}], 1]
                    where disLGX2 uses Conxm[{x,y,z}] := {-x, y, z}
        """
        if step <= self.aL:
            num_r = self.nLr + 1
        else:
            num_r = self.aL + self.lL + 1
        
        num_z = self.HL + 1
        num_theta_half = self.nL_theta + 1
        NLr_theta = num_r * num_theta_half
        
        # Reshape data: [num_z layers, each with (nLθ+1) * (nLr+1) nodes]
        # Then further reshape each layer to [nLθ+1 angles, nLr+1 radii]
        data_extended = []
        
        for iz in range(num_z):
            # Extract one z-layer
            start_idx = iz * NLr_theta
            end_idx = (iz + 1) * NLr_theta
            layer_data = data_half[start_idx:end_idx]
            
            # Reshape to [num_theta_half, num_r, 3]
            layer_reshaped = layer_data.reshape(num_theta_half, num_r, 3)
            
            # Original half: θ = 0° to 90° (indices 0 to nLθ)
            # Mirror half: θ = 90° to 180° (mirror of indices nLθ-1 down to 1)
            # Note: 0° (index 0) and 90° (index nLθ) are not duplicated
            
            # Add original half
            for itheta in range(num_theta_half):
                for ir in range(num_r):
                    data_extended.append(layer_reshaped[itheta, ir])
            
            # Add mirrored half
            # Mathematica: Range[nLθ] with formula (nLθ - #2) where #2 = 1~30
            # gives theta indices: 29, 28, ..., 1, 0 (30 values)
            # In 0-indexed Python: range(29, -1, -1) = 29 down to 0
            for itheta in range(self.nL_theta - 1, -1, -1):  # nLθ-1 down to 0
                for ir in range(num_r):
                    disp = layer_reshaped[itheta, ir].copy()
                    disp[0] = -disp[0]  # Conxm: {-x, y, z}
                    data_extended.append(disp)
        
        return np.array(data_extended)
    
    def load_and_extend_data(self, step):
        """Load data from file and extend to full circle"""
        log_dir = self.result_root / f"step{step:05d}" / "log"
        
        # Load half-circle data (0-90°)
        disLG = np.loadtxt(log_dir / "u_gl.l.dat", skiprows=1)[:, 1:]
        velLG = np.loadtxt(log_dir / "v_gl.l.dat", skiprows=1)[:, 1:]
        acceLG = np.loadtxt(log_dir / "a_gl.l.dat", skiprows=1)[:, 1:]
        
        # Extend to full circle (0-180°)
        disLGX = self.extend_data_to_full_circle(disLG, step)
        velLGX = self.extend_data_to_full_circle(velLG, step)
        acceLGX = self.extend_data_to_full_circle(acceLG, step)
        
        return disLGX, velLGX, acceLGX
    
    def calJ(self, step):
        """Main J-integral calculation for one step"""
        print(f"Step {step}: Calculating J-integral...")
        
        # Generate mesh (nodeLX, elemLX)
        nodeLX = self.generate_mesh(step)
        elemLX = self.generate_elements(step)
        nnmLX = len(nodeLX)
        nemLX = len(elemLX)
        
        # Debug: save mesh data to files
        # debug_dir = Path(f"results/step{step:05d}/debug")
        # debug_dir.mkdir(parents=True, exist_ok=True)
        # 
        # # Save nodeLX (node coordinates)
        # np.savetxt(debug_dir / "nodeLX.dat", nodeLX, 
        #        header=f"Node coordinates for step {step}\nNode_ID X Y Z",
        #        fmt='%.15e')
        
        # Save elemLX (element connectivity, 1-based indexing)
        # np.savetxt(debug_dir / "elemLX.dat", elemLX, 
        #        header=f"Element connectivity for step {step}\nElem_ID N1 N2 N3 N4 N5 N6 N7 N8",
        #        fmt='%d')
        # 
        # print(f"  Debug: Mesh data saved to {debug_dir}")
        
        print(f"  Generated mesh: {nnmLX} nodes, {nemLX} elements")
        
        # Load and extend data
        disLGX, velLGX, acceLGX = self.load_and_extend_data(step)
        
        # Transform to cylindrical coordinates
        r_theta_y = np.zeros((nnmLX, 3))
        for i in range(nnmLX):
            x, y, z = nodeLX[i]
            r = np.sqrt(x**2 + y**2)
            if x == 0.0:
                theta = np.pi / 2
            else:
                theta = np.arctan2(y, x)
            r_theta_y[i] = [r, theta, z]
        
        # Debug: save extended data to files
        # np.savetxt(debug_dir / "disLGX.dat", disLGX,
        #        header=f"Extended displacement data for step {step}\nNode_ID Ux Uy Uz",
        #        fmt='%.15e')
        # np.savetxt(debug_dir / "velLGX.dat", velLGX,
        #        header=f"Extended velocity data for step {step}\nNode_ID Vx Vy Vz",
        #        fmt='%.15e')
        # np.savetxt(debug_dir / "acceLGX.dat", acceLGX,
        #        header=f"Extended acceleration data for step {step}\nNode_ID Ax Ay Az",
        #        fmt='%.15e')
        # np.savetxt(debug_dir / "r_theta_y.dat", r_theta_y,
        #        header=f"Cylindrical coordinates for step {step}\nNode_ID r theta y",
        #        fmt='%.15e')
        # 
        # print(f"  Debug: Extended data saved to {debug_dir}")
        
        # front0: nodal set along crack front
        if step <= self.aL:
            front0 = [step + 1 + i * (self.nLr + 1) for i in range(2 * self.nL_theta + 1)]
        else:
            front0 = [self.aL + 1 + i * (self.nLr + 1) for i in range(2 * self.nL_theta + 1)]
        
        # front: subset from middle to end (48° to 90° for odd case, 45° to 90° for even case)
        # For odd number of angles (nL_theta+1 is odd), start from next angle to avoid duplication
        start_idx = (self.nL_theta + 1) // 2
        front = front0[start_idx:self.nL_theta + 1]
        
        # print(f"  front0: {len(front0)} nodes, front: {len(front)} nodes")
        # print(f"  First front node (1-based): {front[0]}")
        # print(f"  Front node indices: {front[:5]}...{front[-3:]}")
        # print(f"  Front angles: {[r_theta_y[n-1, 1]*180/np.pi for n in [front[0], front[-1]]]}")
        
        # Calculate J for each crack front node
        Jintsd0 = []
        for nL in front:
            result = self.getJ(nL, step, nodeLX, elemLX, disLGX, velLGX, acceLGX, r_theta_y, front0)
            Jintsd0.append(result)
        
        # Mirror results from the computed 45-90 half to 0-90.
        angles = [row[0] for row in Jintsd0]
        J_values = [row[1] for row in Jintsd0]
        angles_mirrored = list(reversed([90.0 - a for a in angles]))
        J_mirrored = list(reversed(J_values))
        if self.nL_theta % 2 == 0:
            angles_full = angles_mirrored[:-1] + angles
            J_full = J_mirrored[:-1] + J_values
        else:
            angles_full = angles_mirrored + angles
            J_full = J_mirrored + J_values
        
        # Store results
        if step == self.stepini:
            self.JAll = [[["step"] + angles_full]]
        self.JAll[0].append([step] + J_full)
        
        # Calculate K_I
        v1 = np.sqrt(((1 - self.nu) * self.EE) / ((1 + self.nu) * (1 - 2 * self.nu) * self.rho))
        v2 = np.sqrt(self.EE / ((1 + self.nu) * 2 * self.rho))
        beta1 = np.sqrt(1 - (self.v / v1)**2)
        beta2 = np.sqrt(1 - (self.v / v2)**2)
        AI = (beta1 * (1 - beta2**2)) / (4 * beta1 * beta2 - (1 + beta2**2)**2)
        
        K_I_values = [np.sqrt((self.EE * J) / ((1 + self.nu) * AI)) if J > 0 else 0.0 for J in J_full]
        
        if step == self.stepini:
            self.KIDAll = [[["step"] + angles_full]]
        self.KIDAll[0].append([step] + K_I_values)
        
        # Print first calculated angle (typically 45°) and 0° for comparison
        # first_angle_idx = 0
        # center_angle_idx = len(angles)
        # print(f"  Step {step}: J[θ={angles[first_angle_idx]:.1f}°] = {J_values[first_angle_idx]:.2e}, "
        #       f"J[θ=0°] = {J_full[center_angle_idx]:.2e}, K_I[0°] = {K_I_values[center_angle_idx]:.2e}")
        
        return Jintsd0
    
    def getJ(self, nL, step, nodeLX, elemLX, disLGX, velLGX, acceLGX, r_theta_y, front0):
        """
        Calculate J-integral for crack front node nL (1-based index)
        """
        nnmLX = len(nodeLX)
        
        # Crack front angle
        theta0 = r_theta_y[nL - 1, 1]  # Convert to 0-based
        
        # Current crack length
        c0 = step * self.hL
        
        # Delta_theta
        Delta_theta = (np.pi * self.hL * step) / (2 * self.nL_theta)
        
        # Transform to crack-tip local coordinates
        xyz = np.zeros_like(nodeLX)
        sty = np.zeros_like(nodeLX)
        for i in range(nnmLX):
            r, theta, z_orig = r_theta_y[i]
            x_local = r * np.cos(theta - theta0) - c0
            y_local = z_orig
            z_local = r * np.sin(theta - theta0)
            xyz[i] = [x_local, y_local, z_local]
            
            s = r - c0
            t = r * (theta - theta0)
            sty[i] = [s, t, z_orig]
        
        # q1i: nodes in elements containing crack front node nL
        elements_with_nL = []
        for ie, elem in enumerate(elemLX):
            if nL in elem:
                elements_with_nL.append(ie)
        
        q1i = set()
        for ie in elements_with_nL:
            for node in elemLX[ie]:
                q1i.add(node)
        
        # Weight function q
        qi = np.zeros(nnmLX)
        
        for i in range(nnmLX):
            # If node i+1 (1-based) is in q1i, set qi = 1 directly
            if (i + 1) in q1i:
                qi[i] = 1.0
                continue
            
            s, t, y_coord = sty[i]
            norm_sy = np.sqrt(s**2 + y_coord**2)
            
            # qR
            if norm_sy <= self.Rj0 * self.hL + 1e-8:
                qR = 1.0
            elif self.Rj0 == self.Rj1:
                qR = 0.0
            elif norm_sy < self.Rj1 * self.hL:
                qR = (norm_sy - self.Rj1 * self.hL) / ((self.Rj0 - self.Rj1) * self.hL)
            else:
                qR = 0.0
            
            # qW
            abs_t_2 = 2 * np.abs(t)
            if abs_t_2 <= self.Wj0 * Delta_theta + 1e-8:
                qW = 1.0
            elif self.Wj0 == self.Wj1:
                qW = 0.0
            elif abs_t_2 < self.Wj1 * Delta_theta:
                qW = (abs_t_2 - self.Wj1 * Delta_theta) / ((self.Wj0 - self.Wj1) * Delta_theta)
            else:
                qW = 0.0
            
            qi[i] = qR * qW
        
        # meas (using 0-based indexing for front0)
        meas = Delta_theta * np.sum([qi[idx - 1] for idx in front0])
        
        # Elements with non-zero q
        qe = np.array([np.sum(qi[elem - 1]) for elem in elemLX])
        nq = np.where(qe > 1e-8)[0]
        
        if len(nq) == 0:
            return [theta0 * 180 / np.pi, 0.0, 0.0, 0.0]
        
        elemq = elemLX[nq]
        
        # Element node coordinates in crack-tip system
        enode = np.array([[xyz[node - 1] for node in elem] for elem in elemq])
        
        # Jacobian: J = Transpose[Outer[#1 . #2 &, Dnn, enode, 1]]
        n_elem = len(nq)
        n_gauss = len(self.Dnn)
        J_matrices = np.zeros((n_elem, n_gauss, 3, 3))
        for ie in range(n_elem):
            for ig in range(n_gauss):
                # Mathematica: J = Transpose[Outer[#1 . #2 &, Dnn, enode, 1]]
                # Python equivalent: J = (Dnn @ enode.T).T = enode @ Dnn.T
                J_matrices[ie, ig] = (enode[ie].T @ self.Dnn[ig].T).T
        
        # detJ with NEGATIVE sign
        detJ = -np.linalg.det(J_matrices)
        
        # bb = J^(-1) · Dnn
        bb = np.zeros((n_elem, n_gauss, 3, 8))
        for ie in range(n_elem):
            for ig in range(n_gauss):
                bb[ie, ig] = np.linalg.inv(J_matrices[ie, ig]) @ self.Dnn[ig]
        
        # Rotate displacements, velocities, and accelerations: uxyz
        dispR = np.zeros_like(disLGX)
        velR = np.zeros_like(velLGX)
        acceR = np.zeros_like(acceLGX)
        for i in range(nnmLX):
            u, v, w = disLGX[i]
            ux = u * np.cos(theta0) + v * np.sin(theta0)
            uz = -u * np.sin(theta0) + v * np.cos(theta0)
            uy = w
            dispR[i] = [ux, uy, uz]

            vu, vv, vw = velLGX[i]
            vx = vu * np.cos(theta0) + vv * np.sin(theta0)
            vz = -vu * np.sin(theta0) + vv * np.cos(theta0)
            vy = vw
            velR[i] = [vx, vy, vz]
            
            au, av, aw = acceLGX[i]
            ax = au * np.cos(theta0) + av * np.sin(theta0)
            az = -au * np.sin(theta0) + av * np.cos(theta0)
            ay = aw
            acceR[i] = [ax, ay, az]
        
        # Element displacements, velocities, and accelerations
        dispqe = np.array([[dispR[node - 1] for node in elem] for elem in elemq])
        velqe = np.array([[velR[node - 1] for node in elem] for elem in elemq])
        acceqe = np.array([[acceR[node - 1] for node in elem] for elem in elemq])
        
        # Calculate strain
        epsilon = np.zeros((n_elem, n_gauss, 6))
        for ie in range(n_elem):
            for ig in range(n_gauss):
                b = bb[ie, ig]
                disp = dispqe[ie]
                epsilon[ie, ig, 0] = b[0] @ disp[:, 0]
                epsilon[ie, ig, 1] = b[1] @ disp[:, 1]
                epsilon[ie, ig, 2] = b[2] @ disp[:, 2]
                epsilon[ie, ig, 3] = b[0] @ disp[:, 1] + b[1] @ disp[:, 0]
                epsilon[ie, ig, 4] = b[0] @ disp[:, 2] + b[2] @ disp[:, 0]
                epsilon[ie, ig, 5] = b[1] @ disp[:, 2] + b[2] @ disp[:, 1]
        
        # Calculate stress
        sigma = np.zeros((n_elem, n_gauss, 6))
        for ie in range(n_elem):
            for ig in range(n_gauss):
                sigma[ie, ig] = self.de @ epsilon[ie, ig]
        
        # Calculate Du
        Du = np.zeros((n_elem, n_gauss, 3))
        for ie in range(n_elem):
            for ig in range(n_gauss):
                Du[ie, ig] = bb[ie, ig, 0] @ dispqe[ie]

        # Calculate Dvel/dx
        Dve = np.zeros((n_elem, n_gauss, 3))
        for ie in range(n_elem):
            for ig in range(n_gauss):
                Dve[ie, ig] = bb[ie, ig, 0] @ velqe[ie]
        
        # Calculate Dq
        qie = np.array([qi[elem - 1] for elem in elemq])
        Dq = np.zeros((n_elem, n_gauss, 3))
        for ie in range(n_elem):
            for ig in range(n_gauss):
                Dq[ie, ig, 0] = bb[ie, ig, 0] @ qie[ie]
                Dq[ie, ig, 1] = bb[ie, ig, 1] @ qie[ie]
                Dq[ie, ig, 2] = bb[ie, ig, 2] @ qie[ie]
        
        # Acceleration at Gauss points
        accee = np.zeros((n_elem, n_gauss, 3))
        velee = np.zeros((n_elem, n_gauss, 3))
        qgpe = np.zeros((n_elem, n_gauss))
        for ie in range(n_elem):
            for ig in range(n_gauss):
                velee[ie, ig] = self.nn[ig] @ velqe[ie]
                accee[ie, ig] = self.nn[ig] @ acceqe[ie]
                qgpe[ie, ig] = self.nn[ig] @ qie[ie]
        
        # J-integral static part
        Jints = 0.0
        for ie in range(n_elem):
            for ig in range(n_gauss):
                eps = epsilon[ie, ig]
                sig = sigma[ie, ig]
                du = Du[ie, ig]
                dq = Dq[ie, ig]
                
                W = 0.5 * np.dot(eps, sig)
                
                term1 = (sig[0]*du[0] + sig[3]*du[1] + sig[4]*du[2] - W) * dq[0]
                term2 = (sig[3]*du[0] + sig[1]*du[1] + sig[5]*du[2]) * dq[1]
                term3 = (sig[4]*du[0] + sig[5]*du[1] + sig[2]*du[2]) * dq[2]
                
                Jints += (term1 + term2 + term3) * detJ[ie, ig] * self.weight[ig]
        
        # J-integral dynamic part
        Jintd = 0.0
        for ie in range(n_elem):
            for ig in range(n_gauss):
                acce = accee[ie, ig]
                vel = velee[ie, ig]
                dve = Dve[ie, ig]
                du = Du[ie, ig]
                dq = Dq[ie, ig]
                qgp = qgpe[ie, ig]
                kinetic_term = -0.5 * self.rho * np.dot(vel, vel) * dq[0]
                inertial_term = self.rho * (np.dot(acce, du) - np.dot(vel, dve)) * qgp
                Jintd += (kinetic_term + inertial_term) * detJ[ie, ig] * self.weight[ig]
        
        # Normalize
        Jints = (2.0 * Jints) / meas
        Jintd = 2.0 * Jintd
        Jint = Jints + Jintd
        

        
        return [theta0 * 180 / np.pi, Jint, Jints, Jintd]
    
    def run(self, output_file=None):
        """
        Run J-integral calculation for all steps.
        
        Args:
            output_file: Base name for output files (default: "J_integral_results.csv")
                        Will create <output_file> for J values and 
                        <output_file_without_ext>_KId.csv for K_I values
        """
        if output_file is None:
            output_file = "J_integral_results.csv"
        
        print(f"Starting J-integral calculation from step {self.stepini} to {self.stepend}")
        
        for step in range(self.stepini, self.stepend + 1):
            self.calJ(step)
        
        # Parse output file path
        output_path = Path(output_file)
        output_dir = output_path.parent
        if output_dir != Path('.'):
            output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate K_I filename
        base_name = output_path.stem
        ext = output_path.suffix
        kid_file = output_dir / f"{base_name}_KId{ext}"
        
        # Save J values
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            for row in self.JAll[0]:
                writer.writerow(row)
        
        # Save K_I values
        with open(kid_file, 'w', newline='') as f:
            writer = csv.writer(f)
            for row in self.KIDAll[0]:
                writer.writerow(row)
        
        print(f"\nResults saved:")
        print(f"  - {output_path}")
        print(f"  - {kid_file}")


class FEMReferenceJIntegralCalculator(JIntegralCalculator):
    """
    Calculate FEM reference DSIF on the same local mesh used by a sweep case.

    The 2D axisymmetric FEM data are interpolated to the 0-90 degree local
    mesh, transformed to Cartesian components, and then mirrored to 0-180
    degrees using the same symmetry as the hS-IGA J-integral calculation.
    """

    def __init__(self, *args, fem_data_folder=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.fem_loader = FEMDataLoader(fem_data_folder=fem_data_folder)
        if not self.fem_loader.load_data(self.v):
            raise FileNotFoundError(f"FEM data for velocity {self.v:g} was not found")

    def load_and_extend_data(self, step):
        nodeL = self.generate_half_mesh(step)

        dis_2d = self.fem_loader.get_displacement_at_step(step)
        vel_2d = self.fem_loader.get_velocity_at_step(step)
        acce_2d = self.fem_loader.get_acceleration_at_step(step)
        if dis_2d is None or vel_2d is None or acce_2d is None:
            raise ValueError(f"FEM u/v/a data are incomplete for step {step}")

        disFEM = self.fem_loader.interpolate_2d_to_3d(dis_2d, nodeL)
        velFEM = self.fem_loader.interpolate_2d_to_3d(vel_2d, nodeL)
        acceFEM = self.fem_loader.interpolate_2d_to_3d(acce_2d, nodeL)

        return (
            self.extend_data_to_full_circle(disFEM, step),
            self.extend_data_to_full_circle(velFEM, step),
            self.extend_data_to_full_circle(acceFEM, step),
        )


if __name__ == "__main__":
    calc = JIntegralCalculator()
    calc.run()
