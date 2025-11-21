import datetime
from const import material_property as mp
from const import const_local_mesh as clm
from const import const_global_mesh as cgm

# Crack geometry (from param.txt)
c = 4.0e-3  # Crack radius [m]
thi = 1.0  # Thickness [m]

# Crack velocity (from param.txt)
V = 1000.0  # Velocity [m/s]

# Step control
step_start = 0  # stepini
step_end = 101  # stepend (-1 means calculate from c/hL)
REstart = 0  # 1: restart from dynamic analysis (set to 0 if stepini==0)

# Calculate stepall
if step_end == -1:
    stepall = round(c / clm.hL)
else:
    stepall = step_end

# Boundary condition type
nbcebc = 1  # 0: nodal force, 1: essential boundary condition

# Domain size for global mesh
lGoutxy = 2 * c
lGoutz = c
WidthG = 8.0e-3
HeightG = 4.0e-3

# Calculate global element size
hG = cgm.mu_G * clm.hL * cgm.rGL

# Calculate number of control points
import numpy as np
nPtsX = int(np.ceil(WidthG / hG))
nPtsY = nPtsX
nPtsZ = int(np.ceil(HeightG / hG))

nofix = 1  # Changing boundary conditions: 0=No, 1=Yes

# # J integral
# Rj0 = 1.5
# Rj1 = 1.51
# Wj0 = 0.5
# Wj1 = 0.51

# HHT method
Alpha = 0.
Beta = 0.25
Gamma = 0.5

# integral point
ngp = 2

nrefLlist = 1  # href=2^nref
inc = 1  # num. of increment (basically 1)

# Nu = 0.35   # Poisson's ratio
# ee = 3.2e9  # Young's modulus [Pa]
# Rho = 1170. # density [m/s^2]

# Rayleigh damping - mass coefficient
Alpha_l = 0.

# Beta_l (Rk) is calculated dynamically in input_generator.py using:
# Beta_l = 2.57 * h * sqrt(rho/E)
# where h = hG (global) or hL (local) depending on islocal flag

OPENMP = 24      # number of openmp thread


TEST_NUMBER_LIST = [4]
TEST_NUMBER = TEST_NUMBER_LIST[0]
INC_LIST = [1]
DIR_NAME_ADD_LIST = ["R0.8_No.1", "R0.8_No.2", "R1.0_No.1", "R1.0_No.2", "R1.2_No.1", "R1.2_No.2"]
DIR_NAME_ADD = DIR_NAME_ADD_LIST[TEST_NUMBER-1]
DIR_NAME_TEST_LIST = ["R0_8__NO1", "R0_8__NO2", "R1_0__NO1", "R1_0__NO2", "R1_2__NO1", "R1_2__NO2"]
DIR_NAME_TEST = "Test"
# DAY = datetime.datetime.now().strftime("%Y-%m-%d")
DAY = "test"

LOCAL_01_LIST = [1]
Local01 = LOCAL_01_LIST[0]
DYNAMIC_01_LIST = [1]
HREF_LIST = [3]
HREF = HREF_LIST[0]
NFR_LIST = [1]

TEST_START = 1
TEST_END = 1
GET_CTOD = 0

UBUNTU_DIR = "PMMA_old"


DOS_OPEN = 2    # dos screen 0: Not open　1:Open and close　2:Keeping open

ABO = 0  # 1:abort just before calculation of solver
USER_NAME = "lab"
REPO_NAME = "sfem_linear"

# STEP_TIME = 1.2915155311516764e-6

OUTPUT_FOLDER = "outputfolder"
CALC_STEP = 1
# JOB_NAME = "step_" + str(STEP_START) + "_exp_" + str(CALC_STEP)
# FOLDER_NAME = JOB_NAME

is_cross = False
isgetctod = False
