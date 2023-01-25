from ece230.crystals.diamondlattice import DiamondLattice

class Silicon(DiamondLattice):
    '''Silicon class
    NUM_ATOMS: integer, # atoms in unit cell
    ATOMIC_MASS: float in g/mol, or in atomic mass units
    '''
    def __init__(self):
        super().__init__(
            lattice_constant=5.43095, 
            atomic_mass=28.0855)


