from abc import abstractmethod
from ece230.crystals.unit_num import UnitNum


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
    AVOGADROS_NUMBER: float = UnitNum(6.02214076e23, "1/mol")

    def __init__(self, atoms=1, effective_atoms=None, atomic_mass=1.0, lattice_constant=1.0):
        self.NUM_ATOMS: int = atoms
        self.EFFECTIVE_ATOMS: float = effective_atoms if effective_atoms else atoms
        self.ATOMIC_MASS: float = UnitNum(atomic_mass, "g/mol")
        self.LATTICE_CONSTANT: float = UnitNum(lattice_constant, "angstroms")
        self.UNIT_CELL_VOLUME: float = UnitNum(1.0, "Å^3")
        self.ATOM_CONCENTRATION: float = UnitNum(self.atom_concentration(), "N/cm^3")
        self.MASS_DENSITY: float = UnitNum(self.mass_density(), "g/cm^3")

    @abstractmethod    
    def nearest_neighbor_dist(self):
        '''Calculate the nearest neighbor distance for a crystal
        '''
        pass

    def atom_concentration(self):
        '''Calculate the atom concentration for a crystal in (N/cm^3)
        '''
        return UnitNum(self.EFFECTIVE_ATOMS / self.UNIT_CELL_VOLUME, "N/cm^3")

    def mass_density(self):
        '''Calculate the mass density for a crystal in (g/cm^3)
        '''
        return UnitNum(self.ATOM_CONCENTRATION * self.ATOMIC_MASS / self.AVOGADROS_NUMBER, "g/cm^3")

