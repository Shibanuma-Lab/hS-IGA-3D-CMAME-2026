"""
Circular Crack in 3D Solid - S-IGA Simulation Framework

This package provides mesh generation and simulation setup for 
circular crack propagation analysis using Smoothed Isogeometric Analysis.
"""

__version__ = "1.0.0"
__author__ = "Converted from Mathematica"

from . import global_mesh
from . import local_mesh
from . import boundary
from . import input_generator
from . import initial
from . import linux_command

__all__ = [
    'global_mesh',
    'local_mesh',
    'boundary',
    'input_generator',
    'initial',
    'linux_command',
]
