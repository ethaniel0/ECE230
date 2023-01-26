from ece230.crystals.unit_cells import DiamondLattice
from ece230.crystals.unit_num import UnitNum

class Silicon(DiamondLattice):
    '''Silicon class
    NUM_ATOMS: integer, # atoms in unit cell
    ATOMIC_MASS: float in g/mol, or in atomic mass units
    '''
    UNIT_CELL_VOLUME = UnitNum(5.43095**3, "Å^3")

    def __init__(self):
        super().__init__(
            lattice_constant=5.43095, 
            atomic_mass=28.0855)



