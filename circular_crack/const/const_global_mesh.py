# Global mesh parameters (from param.txt)

# B-spline degrees
p = 2  # Degree in X direction
q = 2  # Degree in Y direction  
r = 2  # Degree in Z direction

# Mesh refinement parameters
nxy = 5  # Number of subdivisions in XY
rxy = [1.2, 2.0]  # Refinement ratios in XY
nGminz = 17  # Minimum number of elements in Z
rz = [1.2, 2.0]  # Refinement ratios in Z

# Global/local element size ratio
rGL = 6  # Ratio of global to local element size (hG/hL)

# Mesh scaling
mu_G = 0.99 ** 0.5  # Global mesh size multiplier