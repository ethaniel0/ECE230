class Crystal:
    '''A collection of basic crystal and unit cell properties
    NUM_ATOMS: integer, # atoms in unit cell
    EFFECTIVE_ATOMS: float, # effective atoms in unit cell, based on geometry (e.g. 1/8 for a corner atom, 1/4 for an edge atom, 1/2 for a face atom, 1 for a body atom)
    ATOMIC_MASS: float in g/mol, or in atomic mass units
    LATTICE_CONSTANT: float in angstroms, length of unit cell edge
    UNIT_CELL_VOLUME: float in angstroms^3, volume of unit cell
    ATOM_CONCENTRATION: float in N/cm^3, concentration of atoms in unit cell
    MASS_DENSITY: float in g/cm^3, mass density of unit cell
    '''

    def __init__(self, atoms=0, effective_atoms=0, atomic_mass=0.0, lattice_constant=0.0, unit_cell_volume=0.0):
        self.NUM_ATOMS: int = atoms
        self.EFFECTIVE_ATOMS: float = effective_atoms
        self.ATOMIC_MASS: float = atomic_mass
        self.LATTICE_CONSTANT: float = lattice_constant
        self.UNIT_CELL_VOLUME: float = unit_cell_volume
        self.ATOM_CONECNTRATION: float = 0.0
        self.MASS_DENSITY: float = 0.0

    def calc_nearest_neighbor_dist():
        '''Calculate the nearest neighbor distance for a crystal
        '''
        pass

    def calc_atom_concentration(self):
        '''Calculate the atom concentration for a crystal in (N/cm^3)
        '''
        self.ATOM_CONECNTRATION = self.EFFECTIVE_ATOMS / self.UNIT_CELL_VOLUME

    def calc_mass_density(self):
        '''Calculate the mass density for a crystal in (g/cm^3)
        '''
        AVOGADROS_NUMBER = 6.02214076e23
        self.MASS_DENSITY = self.ATOM_CONECNTRATION * self.ATOMIC_MASS / AVOGADROS_NUMBER