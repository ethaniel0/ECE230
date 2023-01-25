from ece230.crystals.crystal import Crystal
import sympy as sp

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

if __name__ == "__main__":
    sc = SimpleCubic()
    print(sc.NUM_ATOMS)

class DiamondLattice(Crystal):
    NUM_ATOMS: int = 18
    EFFECTIVE_ATOMS: float = 8

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



    def calc_center_to_center_dist(self):
        '''Calculate the center to center distance for a diamond lattice
        '''
        # Calculate the center to center distance by making a triangle
        #Uses corner atom and bonding atom
        base = sp.sqrt(2) * self.LATTICE_CONSTANT / 4
        height = self.LATTICE_CONSTANT / 4
        return sp.sqrt(base**2 + height**2)
