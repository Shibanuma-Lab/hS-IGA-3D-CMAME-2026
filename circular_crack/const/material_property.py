# Material properties (from param.txt)

SigmaInfinity = 1.0e11  # Applied stress [Pa]
Nu = 0.3  # Poisson's ratio
EE = 2.06e11  # Young's modulus [Pa]
SigmaY0 = 400.0e6  # Yield stress [Pa] (not used)
Rho = 7800.0  # Density [kg/m^3]

# Static mode values (dimensionless/normalized)
SigmaInfinity_static = 1.0  # Normalized applied stress (dimensionless)
Nu_static = 0.3  # Poisson's ratio
EE_static = 100.0  # Normalized Young's modulus (dimensionless)

def update_for_static_mode():
    """Update material properties for static mode (dimensionless system)"""
    global SigmaInfinity, Nu, EE
    SigmaInfinity = SigmaInfinity_static
    Nu = Nu_static
    EE = EE_static

