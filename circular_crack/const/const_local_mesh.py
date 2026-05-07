# Local mesh parameters (from param.txt)

# Local mesh parameters (from param.txt)
# Note: Values will be updated if static_mode is enabled in simulation_params

# Default values (dynamic mode)
hL = 0.04 * 10**-3  # Length of local element [m]
aL = 16
lL = 12
HL = 15
d_theta = 2.604577870356  # Angular resolution [degrees]

# Static mode values (will be set by update_for_static_mode)
# Note: Static mode uses dimensionless (normalized) units
#       All values are pure numbers without physical units
hL_static = 0.011363636363636364  # Normalized element size (dimensionless)
aL_static = 22
lL_static = 22
HL_static = 22

def update_for_static_mode():
    """Update parameters for static mode"""
    global hL, aL, lL, HL
    hL = hL_static
    aL = aL_static
    lL = lL_static
    HL = HL_static
