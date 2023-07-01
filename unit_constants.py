from sympy.physics import units
from math import pi

angstroms = 1e-8 * units.cm
cal = 0.239006 * units.joules

a_si: units.Quantity = 5.43086 * angstroms
'''
lattice constant of silicon

5.43086 Å
'''

N_si: units.Quantity = 5e22 * units.cm**-3
'''
# silicon atoms per cubic cm

5e22 cm^-3
'''

e_thermal_vel: units.Quantity = 1e7 * units.cm / units.seconds
'''
root mean square velocity of elections in silicon at 300K

1e7 cm/s
'''

q: units.Quantity = 1.602177e-19 * units.coulomb
'''
electron charge
1.602177e-19 C
'''

Kb: units.Quantity = 1.38058e-23 * units.joule / units.kelvin
'''
Boltzmann's Constant
1.38058e-23 J/K
'''

KbToq: units.Quantity = 0.025852 * units.volts
'''
KbT/q at 300K

0.025852 V
'''

KbT: units.Quantity = KbToq * q
'''
KbT at 300K

4.1419479804e-21 coulomb
'''


h: units.Quantity = 4.135669e-15 * units.electronvolts * units.seconds
'''
Planck's Constant
4.135669e-15 eV*s
'''

hbar: units.Quantity = h / (2 * pi)
'''
Reduced Planck's Constant

6.58212164341916e-16 eV*s
'''

m0: units.Quantity = 9.109389e-28 * units.grams
'''
Electron Rest Mass 
9.109389e-28 grams
'''

mn_300k: units.Quantity = 0.2588 * m0
'''
electron conductivity effective mass at 300K

also can be the electron effective mass
'''

mp_300k: units.Quantity = 0.38158 * m0
'''
hole conductivity effective mass at 300K

also can be the hole effective mass
'''

m_ndos_300K: units.Quantity = 1.09*m0
'''
electron density-of-states effective mass 

9.92923401e-28 grams
'''

m_pdos_300K: units.Quantity = 1.15*m0
'''
hole density-of-states effective mass 

1.047579735e-27 grams
'''

m_ndos_400K: units.Quantity = 1.11*m0
'''
electron density-of-states effective mass 

1.011142179e-27 grams
'''

m_pdos_400K: units.Quantity = 1.23*m0
'''
hole density-of-states effective mass 

1.120454847e-27 grams
'''

c: units.Quantity = 2.99792458e10 * (units.cm / units.seconds)
'''
Speed of light

c = 2.99792458e10 cm/s
'''

eps0: units.Quantity = 8.854187e-14 * (units.farads / units.cm)
'''
Vacuum Permittivity

8.854187e-14 F/cm
'''

mu0: units.Quantity = 1.256637e-8 * (units.henrys / units.cm)
'''
Vacuum Permeability

1.256637e-8 H/cm
'''

Navo: units.Quantity = 6.022136e23 * (1 / units.moles)
'''
Avogadro's Number

6.022136e23 mol^-1
'''

R: units.Quantity = 1.98719 * cal/(units.mole*units.kelvin)
'''
Gas Constant

1.98719 cal/(mol*K)
'''

Vm: units.Quantity = 22.41410 * units.liters/units.moles
'''
Molar Volume

22.41410 lit/mol
'''

ab: units.Quantity = 0.529177 * angstroms
'''
Bohr's Radius

0.529177 Å
'''

ni: units.Quantity = 1.07e10 * units.cm**(-3)
'''
Si Intrinsic Carrier Conc. at 300 K
'''

eps_si: units.Quantity = 11.7 * eps0
'''
Silicon Permittivity

1.035939879e-12 F/cm
'''

eps_ox: units.Quantity = 3.9 * eps0
'''
Silicon Dioxide Permittivity

3.45313293e-13 F/cm
'''

Nc: units.Quantity = 2.86e19 * units.cm**-3
'''
conduction-band effective density of states at 300 K

2.86e19 states/cm^-3
'''

Nv: units.Quantity = 3.1e19 * units.cm**-3
'''
valence-band effective density of states at 300 K

3.1e19 states/cm^-3
'''

Nsi: units.Quantity = 5e22 * units.cm**-3
'''
Silicon atomic density

5e22 atoms/cm^3
'''

Eg_si: units.Quantity = 1.1242 * units.electronvolts
'''
Silicon band gap energy at 300K

1.1242 eV
'''

eg_si0: units.Quantity = 1.17 * units.electronvolts
'''
Silicon band gap energy at 0K

1.17 eV
'''

electron_sat_vel_300K: units.Quantity = 1.066e7 * units.cm/units.seconds
'''
saturation velocity of electrons in silicon at 300K

1.066e7 * cm/s
'''

hole_sat_vel_300K: units.Quantity = 8.29e6 * units.cm/units.seconds
'''
saturation velocity of holes in silicon at 300K

8.29e6 * cm/s
'''