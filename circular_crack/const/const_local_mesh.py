# Local mesh parameters (from param.txt)

# Local mesh parameters (from param.txt)
# Note: Values will be updated if static_mode is enabled in simulation_params

# Default values (dynamic mode)
hL = 0.04 * 10**-3  # Length of local element [m]
aL = 20  # Number of local elements on crack surface
lL = 15  # Number of local elements on ligament surface
HL = 18  # Number of local elements in thickness direction
d_theta = 3  # Angular resolution [degrees]

# Static mode values (will be set by update_for_static_mode)
# Note: Static mode uses dimensionless (normalized) units
#       All values are pure numbers without physical units
hL_static = 1.0 / 48.0  # Normalized element size (dimensionless)
aL_static = 12
lL_static = 12
HL_static = 12

def update_for_static_mode():
    """Update parameters for static mode"""
    global hL, aL, lL, HL
    hL = hL_static
    aL = aL_static
    lL = lL_static
    HL = HL_static
