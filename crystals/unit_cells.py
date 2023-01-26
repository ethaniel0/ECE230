from ece230.crystals.crystal import Crystal
from ece230.crystals.unit_num import UnitNum
import sympy as sp
import numpy as np

class SimpleCubic(Crystal):
    NUM_ATOMS: int = 8
    EFFECTIVE_ATOMS: float = 1
    
    def __init__(self, lattice_constant=1.0) -> None:
        '''A simple cubic crystal, inheriting from Crystal
        NUM_ATOMS: 8 atoms in unit cell
        EFFECTIVE_ATOMS: 1 effective atom in unit cell, based on geometry (e.g. 1/8 for a corner atom, 1/4 for an edge atom, 1/2 for a face atom, 1 for a body atom)
        ATOMIC_MASS: float in g/mol, or in atomic mass units
        LATTICE_CONSTANT: float in angstroms, length of unit cell edge
        UNIT_CELL_VOLUME: float in angstroms^3, volume of unit cell
        '''
        
        self.__doc__ = Crystal.__doc__
        super().__init__(
            atoms=SimpleCubic.NUM_ATOMS, 
            effective_atoms=SimpleCubic.EFFECTIVE_ATOMS, 
            lattice_constant=lattice_constant)
        
        self.UNIT_CELL_VOLUME = UnitNum(lattice_constant**3, "Å^3")
    
    def nearest_neighbor_dist(self) -> sp.Rational:
        '''Calculate the nearest neighbor distance for a simple cubic crystal
        '''
        return self.LATTICE_CONSTANT

    


class DiamondLattice(Crystal):
    NUM_ATOMS: int = 18
    EFFECTIVE_ATOMS: float = 8
    BOND_ANGLE: float = UnitNum(sp.acos(-1/3), 'rad')

    def __init__(self, lattice_constant=1.0, atomic_mass=1) -> None:
        '''A diamond lattice crystal, inheriting from Crystal
        NUM_ATOMS: 18 atoms in unit cell
        EFFECTIVE_ATOMS: 8 effective atom in unit cell, based on geometry (e.g. 1/8 for a corner atom, 1/4 for an edge atom, 1/2 for a face atom, 1 for a body atom)
        ATOMIC_MASS: float in g/mol, or in atomic mass units
        LATTICE_CONSTANT: float in angstroms, length of unit cell edge
        UNIT_CELL_VOLUME: float in angstroms^3, volume of unit cell
        '''
        self.__doc__ = Crystal.__doc__
        super().__init__(
            atoms=DiamondLattice.NUM_ATOMS, 
            effective_atoms=DiamondLattice.EFFECTIVE_ATOMS, 
            lattice_constant=lattice_constant,
            atomic_mass=atomic_mass)
        
        self.UNIT_CELL_VOLUME = self.LATTICE_CONSTANT**3

    def center_to_center_dist(self) -> sp.Rational:
        '''Calculate the center to center distance between atoms for a diamond lattice
        '''
        # Calculate the center to center distance by making a triangle
        #Uses corner atom and bonding atom
        base = sp.sqrt(2) * sp.Rational(self.LATTICE_CONSTANT, 4)
        height = sp.Rational(self.LATTICE_CONSTANT, 4)
        return UnitNum(sp.sqrt(base**2 + height**2), 'Å')
