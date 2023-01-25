import sympy as sp

def calc_miller_inds(x: float, y: float, z: float) -> tuple[int, int, int]:
    """Calculate the Miller indices for a crystal given the lattice parameters.
    Args:
        x: The x intercept.
        y: The y intercept.
        z: The z intercept.
    Returns:
        A tuple of the Miller indices.

    Process:
        1. Indices h, k, l = 1/x, 1/y, and 1/z respectively
        2. Multiply each index by the smallest number such that all indices are integers
    """
    x, y, z = sp.Rational(x), sp.Rational(y), sp.Rational(z)
    h , k, l = 1/x, 1/y, 1/z

    min_val = min(h, k, l)
    h, k, l = h/min_val, k/min_val, l/min_val

    # get the lcm of the denominators of the miller indices, multiply by that
    lcm_num = sp.lcm([h.q, k.q, l.q])
    h, k, l = lcm_num*h, lcm_num*k, lcm_num*l

    return (int(h), int(k), int(l))

def calc_miller_inds_decimal(x: float, y: float, z: float) -> tuple[float, float, float]:
    """Calculate the Miller indices for a crystal given the lattice parameters.
    Args:
        x: The x lattice parameter.
        y: The y lattice parameter.
        z: The z lattice parameter.
    Returns:
        A tuple of the Miller indices.
    
    Process:
        1. Indices h, k, l = 1/x, 1/y, and 1/z respectively
        2. Multiply each index by the smallest number such that all indices are integers
    """
    h , k, l = 1/x, 1/y, 1/z

    min_val = min(h, k, l)
    h, k, l = h/min_val, k/min_val, l/min_val

    return (h, k, l)

if __name__ == "__main__":
    print(calc_miller_inds(9.66, 19.32, 14.49))
    print(calc_miller_inds_decimal(9.66, 19.32, 14.49))