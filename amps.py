from solvers.mosfet import Mosfet
from sympy import symbols, solve

def ssCommonQPoint(Vdd: float, Rd: float, Rs: float, mf: Mosfet):
    '''
    Calculates the Q-point for a common single-stage amplifier.
    '''
    Vds = symbols("Vds")
    
    id1 = (Vdd - Vds) / (Rd + Rs)
    id2 = mf.drain_current_lin_kn(mf.Kn, mf.Vgs, mf.Vtn, Vds)
    
    VdsQ = solve(id1 - id2, Vds)[0]
    
    print(VdsQ)
    
    if VdsQ > mf.Vdsat:
        id2 = mf.drain_current_sat_kn(mf.Kn, mf.Vgs, mf.Vtn, mf.ch_len_modulation, Vds, mf.Vdsat)
        VdsQ = solve(id1 - id2, Vds)[0]
    
    IdsQ = id1.subs(Vds, VdsQ)
    
    return IdsQ, VdsQ