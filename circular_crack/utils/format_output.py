"""
Output formatting utilities to match Mathematica format exactly
"""

def format_real(x):
    """
    Format a real number to match Mathematica's output format
    
    - 15 significant digits total (14 decimal places + 1 before decimal)
    - Scientific notation with uppercase 'E'
    - Always in exponential form: 1.23456789012345E+00
    - Zero is formatted as 0.000000000000000E+00
    
    Args:
        x: float number to format
    
    Returns:
        Formatted string matching Mathematica format
    """
    # Format with 14 decimal places (total 15 significant digits)
    # Python's .14E gives: 1.23456789012345E+00 (1 + 14 = 15 digits)
    formatted = f"{x:.15E}"
    return formatted


def format_node_line(node_id, x, y, z):
    """
    Format a node line for node.g.dat or node.l.dat
    
    Args:
        node_id: Node ID (integer)
        x, y, z: Coordinates (float)
    
    Returns:
        Formatted string with tab separators
    """
    return f"{node_id}\t{format_real(x)}\t{format_real(y)}\t{format_real(z)}"


def format_bc_line(node_id, dof, value):
    """
    Format a boundary condition line for bc.g.dat or bc.l.dat
    
    Args:
        node_id: Node ID (integer)
        dof: Degree of freedom (integer, 1-3)
        value: BC value (float)
    
    Returns:
        Formatted string with tab separators
    """
    return f"{node_id}\t{dof}\t{format_real(value)}"


def format_weight_line(node_id, weight):
    """
    Format a weight line for weights.g.dat
    
    Args:
        node_id: Node ID (integer)
        weight: Weight value (float)
    
    Returns:
        Formatted string with tab separators
    """
    return f"{node_id}\t{format_real(weight)}"


def format_init_line(node_id, dof, value):
    """
    Format an initial condition line for delta_u.dat, v_init.dat, etc.
    
    Args:
        node_id: Node ID (integer)
        dof: Degree of freedom (integer, 1-3)
        value: Initial value (float)
    
    Returns:
        Formatted string with tab separators
    """
    return f"{node_id}\t{dof}\t{format_real(value)}"
