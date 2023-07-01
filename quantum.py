from solvers import constants as cs
from sympy import exp, symbols, pi, sqrt

def fermiDirac(E=symbols("E"), EF=symbols("EF"), KbT=symbols("KbT"), numerical=False):
    '''
    Fermi-Dirac Equation, calculates distribution of particles among quantum states.

    1/(1 + exp((E-EF)/KbT))
    '''
    if numerical:
        KbT = cs.KbToq
    return 1/(1 + exp((E-EF)/KbT))

def boltzmann(E=symbols("E"), EF=symbols("EF"), KbT=symbols('KbT'), numerical=False):
    '''
    Boltzman approximation of the Fermi-Dirac equation
    '''
    if numerical:
        KbT = cs.KbToq
    return exp(-(E-EF)/KbT)

def gc(E=symbols("E"), Ec=symbols("Ec"), mndos=symbols("mndos"), h=symbols("h"), numerical=False):
    '''
    density of quantum states in the conduction band per unit volume
    
    (1/2π^2) * (2*mndos/hbar^2)^(3/2) * sqrt(E - Ec)
    '''
    if numerical:
        mndos = cs.m_ndos_300K
        h = cs.hbar
        
    p1 = 1/(2*pi**2)
    p2 = (2*mndos/(h**2))**(3/2)
    p3 = sqrt(E - Ec)
    return p1 * p2 * p3

def gv(E=symbols("E"), Ev=symbols("Ev"), mpdos=symbols("mpdos"), numerical=False):
    '''
    density of quantum states in the valence band per unit volume
    
    (1/2π^2) * (2*mndos/hbar^2)^(3/2) * sqrt(Ev - E)
    '''
    if numerical:
        mpdos = cs.m_pdos_300K
        
    p1 = 1/(2*pi**2)
    p2 = (2*mpdos/(cs.hbar**2))**(3/2)
    p3 = sqrt(Ev - E)
    return p1 * p2 * p3
    
    
if __name__ == "__main__":
    print(fermiDirac())