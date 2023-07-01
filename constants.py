from sympy.physics import units
from math import pi

a_si: float = 5.43086
'''
lattice constant of silicon

5.43086 Å
'''

N_si: float = 5e22
'''
# silicon atoms per cubic cm

5e22 cm^-3
'''

e_thermal_vel: float = 1e7
'''
root mean square velocity of elections in silicon at 300K

1e7 cm/s
'''

q: float = 1.602177e-19
'''
electron charge
1.602177e-19 C
'''

Kb: float = 1.38058e-23
'''
Boltzmann's Constant
1.38058e-23 J/K
'''

h: float = 4.135669e-15
'''
Planck's Constant
4.135669e-15 eV*s
'''

hbar: float = h / (2 * pi)
'''
Reduced Planck's Constant

6.58212164341916e-16 eV*s
'''

m0: float = 9.109389e-28
'''
Electron Rest Mass 
9.109389e-28 grams
'''

mn_300k: float = 0.2588 * m0
'''
electron conductivity effective mass at 300K

also can be the electron effective mass
'''

mp_300k: float = 0.38158 * m0
'''
hole conductivity effective mass at 300K

also can be the hole effective mass
'''


m_ndos_300K: float = 1.09*m0
'''
electron density-of-states effective mass 

9.92923401e-28 grams
'''

m_pdos_300K: float = 1.15*m0
'''
hole density-of-states effective mass 

1.047579735e-27 grams
'''

m_ndos_400K: float = 1.11*m0
'''
electron density-of-states effective mass 

1.011142179e-27 grams
'''

m_pdos_400K: float = 1.23*m0
'''
hole density-of-states effective mass 

1.120454847e-27 grams
'''

c: float = 2.99792458e10
'''
Speed of light

c = 2.99792458e10 cm/s
'''

eps0: float = 8.854187e-14
'''
Vacuum Permittivity

8.854187e-14 F/cm
'''

mu0: float = 1.256637e-8
'''
Vacuum Permeability

1.256637e-8 H/cm
'''

Navo: float = 6.022136e23
'''
Avogadro's Number

6.022136e23 mol^-1
'''

R: float = 1.98719
'''
Gas Constant

1.98719 cal/(mol*K)
'''

Vm: float = 22.41410
'''
Molar Volume

22.41410 lit/mol
'''

ab: float = 0.529177
'''
Bohr's Radius

0.529177 Å
'''

KbT: float = 0.025852 * q
'''
KbT at 300K

4.1419479804e-21 coulomb
'''

KbToq: float = 0.025852
'''
KbT at 300K

0.025852 V
'''

ni: float = 1.07e10
'''
Si Intrinsic Carrier Conc. at 300 K
'''

eps_si: float = 11.7 * eps0
'''
Silicon Permittivity

1.035939879e-12 F/cm
'''

eps_ox: float = 3.9 * eps0
'''
Silicon Dioxide Permittivity

3.45313293e-13 F/cm
'''

Nc: float = 2.86e19
'''
conduction-band effective density of states at 300 K

2.86e19 states/cm^-3
'''

Nv: float = 3.1e19
'''
valence-band effective density of states at 300 K

3.1e19 states/cm^-3
'''

Nsi: float = 5e22
'''
Silicon atomic density

5e22 atoms/cm^3
'''

Eg_si: float = 1.1242
'''
Silicon band gap energy at 300K

1.1242 eV
'''

eg_si0: float = 1.17
'''
Silicon band gap energy at 0K

1.17 eV
'''

electron_sat_vel_300K: float = 1.066e7
'''
saturation velocity of electrons in silicon at 300K

1.066e7 * cm/s
'''

hole_sat_vel_300K: float = 8.29e6
'''
saturation velocity of holes in silicon at 300K

8.29e6 cm/s
'''