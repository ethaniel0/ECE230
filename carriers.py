from solvers import constants as cs
from solvers.quantum import fermiDirac, gc, boltzmann
from sympy import symbols, integrate, sqrt

def n_thermal_eq_quantum(Ec=symbols('Ec'), EF=symbols('EF'), Evac=symbols('Evac'), numerical=True):
    '''
    calculates electron concentration n◦ at thermal equilibrium in a silicon region using Boltzmann approximation
    
    - Ec = conduction band energy
    - EF = Fermi energy
    - Evac = vacuum energy
    - numerical = replace known constants with numerical values
    '''
    E = symbols('E')
    expr = gc(E, Ec, numerical=numerical) * boltzmann(E, EF, numerical=numerical)
    return integrate(expr, (E, Ec, Evac))
    
def n_thermal_eq_nc(Nc=symbols('Nc'), Ec=symbols("Ec"), EF=symbols('EF'), numerical=True):
    '''
    calculates electron concentration n◦ at thermal equilibrium in a silicon region using classical approximation
    
    - Nc = conduction band effective density of states
    - Ec = conduction band energy
    - EF = Fermi energy
    - numerical = replace known constants with numerical values
    '''
    if numerical:
        Nc = cs.Nc
    return Nc * boltzmann(Ec, EF, numerical=numerical)

def p_thermal_eq_nc(Nv=symbols('Nv'), EF=symbols("EF"), Ev=symbols('Ev'), numerical=True):
    '''
    calculates electron concentration n◦ at thermal equilibrium in a silicon region using classical approximation
    
    - Nv = valence band effective density of states
    - EF = Fermi energy
    - Ev = valence band edge
    - numerical = replace known constants with numerical values
    '''
    if numerical:
        Nv = cs.Nv
    return Nv * boltzmann(EF, Ev, numerical=numerical)

def rec_life_p(Na):
    '''
    minority recombination lifetime of electrons
    
    (τ_rec,P or τ_n,P)
    
    - Na = acceptor concentration
    '''
    denom = 3.45e-12 * Na + 9.5e-32 * (Na**2)
    return 1/denom

def rec_life_n(Nd):
    '''
    minority recombination lifetime of holes
    
    (τ_rec_N or τ_p,N)
    
    - Nd = donor concentration
    '''
    denom = 7.8e-13 * Nd + 1.8e-31 * (Nd**2)
    return 1/denom

def gen_life_p(Na):
    '''
    minority generation lifetime of holes
    
    (τ_gen,P)
    
    - Na = acceptor concentration
    '''
    return 75*rec_life_p(Na)

def gen_life_n(Nd):
    '''
    minority generation lifetime of electrons 
    
    (τ_gen,N or τ_p,N)
    
    - Nd = donor concentration
    '''
    return 75*rec_life_n(Nd)

def mu_p_maj(Na):
    '''
    hole mobility (µ_p) in p-type region
    
    - Na = acceptor concentration
    '''
    return 49.7 + 418.3/(1 + (Na/(1.6e17))**0.7)

def mu_p_min(Nd):
    '''
    hole mobility (µ_p) in n-type region
    
    - Nd = donor concentration
    '''
    return 130 + 370/(1 + (Nd/(8e17))**1.25) 

def mu_n_maj(Nd):
    '''
    electron mobility (µ_n) in n-type region
    
    - Nd = donor concentration
    '''
    return 92 + 1268/(1 + (Nd/(1.3e17))**0.91)

def mu_n_min(Na):
    '''
    calculates electron mobility (µ_n) in p-type region
    
    Na = acceptor concentration
    '''
    return 232 + 1180/(1 + (Na/(8e16))**0.9)

def diffusivity(mu):
    '''
    calculates diffusivity (D) based on mobility (mu)
    
    - mu = µ = mobility
    '''
    return cs.KbToq * mu

def diffusion_length(D, tau):
    return sqrt(D*tau)

def v_drift_n(mu, E):
    '''
    low-field bulk electron drift velocity
    
    - mu = µ_n = electron mobility
    - E = electric field strength
    '''
    return -mu*E

def v_drift_p(mu, E):
    '''
    low-field bulk hole drift velocity
    
    - mu = µ_p = hole mobility
    - E = electric field strength
    '''
    return mu*E
    
def J_drift_n(mu, n, E):
    '''
    low-field bulk hole drift current density
    
    - mu = µ_p = hole mobility
    - n = electron concentration
    - E = electric field strength
    '''
    return cs.q * mu * n * E

def J_drift_p(mu, p, E):
    '''
    low-field bulk hole drift current density
    
    - mu = µ_p = hole mobility
    - p = hole concentration
    - E = electric field strength
    '''
    
    return cs.q * mu * p * E
    
def hole_diffusion_flux(D, dp):
    '''
    hole diffusion flux
    
    - D = hole diffusion coefficient
    - dp = hole concentration gradient
    '''
    return -D * dp

def electron_diffusion_flux(D, dn):
    '''
    electron diffusion flux
    
    - D = electron diffusion coefficient
    - dn = electron concentration gradient
    '''
    return D * dn

def J_diffusion_p(D, dp):
    '''
    hole diffusion current density
    
    - D = hole diffusion coefficient
    - dp = hole concentration gradient
    '''
    return cs.q * hole_diffusion_flux(D, dp)

def J_diffusion_n(D, dn):
    '''
    electron diffusion current density
    
    - D = electron diffusion coefficient
    - dn = electron concentration gradient
    '''
    return cs.q * electron_diffusion_flux(D, dn)

def J_n(mu, n, E, D, dn):
    '''
    total electron current density
    
    - mu = µ_n = electron mobility
    - n = electron concentration
    - E = electric field strength
    - D = electron diffusion coefficient
    - dn = electron concentration gradient
    '''
    return J_drift_n(mu, n, E) + J_diffusion_n(D, dn)

def J_p(mu, p, E, D, dp):
    '''
    total hole current density
    
    - mu = µ_p = hole mobility
    - p = hole concentration
    - E = electric field strength
    - D = hole diffusion coefficient
    - dp = hole concentration gradient
    '''
    return J_drift_p(mu, p, E) + J_diffusion_p(D, dp)

def J_n_from_drift_and_diffusion(drift_current, diffusion_current):
    '''
    total electron current density
    
    - drift_current = J_drift_n
    - diffusion_current = J_diffusion_n
    '''
    return drift_current + diffusion_current

def J_p_from_drift_and_diffusion(drift_current, diffusion_current):
    '''
    total hole current density
    
    - drift_current = J_drift_p
    - diffusion_current = J_diffusion_p
    '''
    return drift_current + diffusion_current

def conductivity(mu_n, mu_p):
    '''
    calculates conductivity (σ) based on electron and hole mobility (mu_n, mu_p)
    
    - mu_n = µ_n = electron mobility
    - mu_p = µ_p = hole mobility
    '''
    return cs.q * (mu_n + mu_p)

def resistivity(mu_n, mu_p):
    '''
    calculates resistivity (ρ) based on electron and hole mobility (mu_n, mu_p)
    
    - mu_n = µ_n = electron mobility
    - mu_p = µ_p = hole mobility
    '''
    return 1 / conductivity(mu_n, mu_p)

def resistance(mu_n, mu_p, L, A):
    '''
    calculates resistance (R) based on electron and hole mobility (mu_n, mu_p), length (L) and area (A)
    
    - mu_n = µ_n = electron mobility
    - mu_p = µ_p = hole mobility
    - L = length
    - A = area
    '''
    return resistivity(mu_n, mu_p) * L / A

def resistance_from_resistivity(p, L, A):
    '''
    calculates resistance (R) based on resistivity (p), length (L) and area (A)
    
    - p  = ρ = resistivity
    - L = length
    - A = area
    '''
    return p * L / A

def debye_length(D, tau):
    '''
    calculates the Debye length (L_D) based on diffusion coefficient (D) and relaxation time (tau)
    
    - D = diffusion coefficient
    - tau = dielectric relaxation time
    '''
    return sqrt(D*tau)