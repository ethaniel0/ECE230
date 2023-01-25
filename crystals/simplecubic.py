from ece230.crystals.crystal import Crystal

class SimpleCubic(Crystal):
   
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
            atoms=8, 
            effective_atoms=1, 
            lattice_constant=lattice_constant)

if __name__ == "__main__":
    sc = SimpleCubic()
    print(sc.NUM_ATOMS)