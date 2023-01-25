class Crystal:
    '''A collection of basic crystal and unit cell properties
    NUM_ATOMS: integer, # atoms in unit cell
    EFFECTIVE_ATOMS: float, # effective atoms in unit cell, based on geometry (e.g. 1/8 for a corner atom, 1/4 for an edge atom, 1/2 for a face atom, 1 for a body atom)
    ATOMIC_MASS: float in g/mol, or in atomic mass units
    '''
    NUM_ATOMS: int = 0
    EFFECTIVE_ATOMS: float = 0
    ATOMIC_MASS: float = 0.0

    def calc_nearest_neighbor_dist():
        '''Calculate the nearest neighbor distance for a crystal
        '''
        pass
