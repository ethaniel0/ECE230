from solvers.constants import KbToq, KbT, ni, q, eps_si
import solvers.constants as cs
from sympy import ln, sqrt, exp, symbols
from solvers.silicon import Silicon

def Vbi(Na, Nd):
    '''Returns the built-in voltage of a diode with Na and Nd acceptor / donors concentrations.'''
    return KbToq * ln(Na*Nd/(ni**2))

def W_depl(Na, Nd, Vpn=0):
    '''Returns the width of the depletion region.'''
    vbi = Vbi(Na, Nd)
    part1 = 2*eps_si / q
    part2 = (Na + Nd) / (Na * Nd)
    return sqrt(part1 * part2 * (vbi - Vpn))

def W_depl_p(Na, Nd, Vpn=0):
    '''Returns the width of the p-side depletion region.'''
    vbi = Vbi(Na, Nd)
    part1 = 2*eps_si / q
    part2 = Nd / (Na * (Na + Nd))
    return sqrt(part1 * part2 * (vbi - Vpn))

def W_depl_n(Na, Nd, Vpn=0):
    '''Returns the width of the n-side depletion region.'''
    vbi = Vbi(Na, Nd)
    part1 = 2*eps_si / q
    part2 = Na / (Nd * (Na + Nd))
    return sqrt(part1 * part2 * (vbi - Vpn))

def E_p_depletion(Na, Wdp, x):
    '''Returns the electric field in the p-side depletion region at x away from center.'''
    return -(q*Na/eps_si) * (x + Wdp)

def E_n_depletion(Nd, Wdn, x):
    '''Returns the electric field in the n-side depletion region at x away from center.'''
    return -(q*Nd/eps_si) * (Wdn - x)

def depl_min_conc_at_edge(p_or_n, Vpn=0):
    '''Returns the minority carrier concentration at the edge of the depletion region'''
    return p_or_n * exp(q*Vpn/KbT)
    
def hole_cur_desnsity_at_n_depl_edge(Nd, DpN, LpN, Vpn=0):
    '''Returns the hole current density at the edge of the n-side depletion region.'''
    return ((q*DpN*(ni**2))/(LpN*Nd)) * (exp(q*Vpn/KbT) - 1)

def electron_cur_density_at_p_depl_edge(Na, DnP, LnP, Vpn=0):
    '''Returns the electron current density at the edge of the p-side depletion region.'''
    return ((q*DnP*(ni**2))/(LnP*Na)) * (exp(q*Vpn/KbT) - 1)

def v_n_depletion_edge(Nd, Wdp, Vpn=0):
    '''Returns the electric potential in the n-side depletion region at x away from center.'''
    return q*Nd/(2*eps_si) * Wdp**2

def v_p_depletion_edge(Na, Wdn, Vpn=0):
    '''Returns the electric potential in the p-side depletion region at x away from center.'''
    return q*Na/(2*eps_si) * Wdn**2

class Diode:
    def __init__(self):
        self.p: Silicon = Silicon()
        self.n: Silicon = Silicon()
        
        self.vbi = symbols('vbi')
        '''Built-in voltage'''
        
        self.v_bias = symbols('v_bias')
        '''Bias voltage'''
        
        self.vbi_n = symbols('vbi_n')
        '''Built-in voltage on the n-side'''
        self.vbi_p = symbols('vbi_p')
        '''Built-in voltage on the p-side'''
        
        self.w_depl = symbols('w_depl')
        '''Width of the depletion region'''
        self.w_depl_p = symbols('w_depl_p')
        '''Width of the p-side depletion region'''
        self.w_depl_n = symbols('w_depl_n')
        '''Width of the n-side depletion region'''
        
        self.depl_edge_p_conc = symbols('p_depl_edge')
        '''Minority carrier (hole) concentration at the edge of the n-side depletion region. units: cm^-3'''
        self.depl_edge_n_conc = symbols('n_depl_edge')
        '''Minority carrier (electron) concentration at the edge of the p-side depletion region. units: cm^-3'''
        
        self.depl_edge_n_J = symbols('J_n_depl')
        '''Electron diffusion current density at the edge of the p-side depletion region. units: A/cm^2'''
        self.depl_edge_p_J = symbols('J_p_depl')
        '''Hole diffusion current density at the edge of the n-side depletion region. units: A/cm^2'''
        
        self.e_bi_max = symbols('ε_x_max')
        '''Maximum electric field in the depletion region'''
        
        self.t_rec_scr = symbols('t_rec_scr')
        '''space-charge-region recombination lifetime'''
        
        self.t_gen_scr = symbols('t_rec_scr')
        '''space-charge-region generation lifetime'''
        
        self.Rmax_rec = symbols('Rmax_rec')
        '''R^net_rec max, electron-hole-pair recombination rate in the depletion region'''
        
        self.J_d_scr = symbols('J_D_scr')
        '''J_d_scr, space-charge-region recombination current density in forward bias'''
        
        self.J_d_scg = symbols('J_D_scg')
        '''J_d_scg, space-charge-generation current density under reverse bias'''
        
        self.J_S_diff = symbols('J_S_diff')
        '''J_S_diff, diffusion saturation current density'''
        
        self.J_S_scr = symbols('J_S_scr')
        '''J_S_scr, space-charge-recombination saturation current density'''
        
        self.J_D = symbols('J_D')
        '''J_D, diode current density'''
        
        self.area = symbols('A')
        
        self.I_D = symbols('I_D')
        
        self.V_breakdown = symbols('V_B')
        '''Breakdown voltage'''
        
        self.Q_pn_dep = symbols('Q_pn_dep')
        '''the depletion charge density'''
        
        self.Q_pn_diff_n = symbols('Q_pn_diff_n')
        '''hole diffusion cahrge density in a long-base n-side quasi-netral region'''
        
        self.Q_pn_diff_p = symbols('Q_pn_diff_p')
        '''electron diffusion cahrge density in a long-base p-side quasi-netral region'''
        
        self.C_pn_dep = 0
        '''small-signal depletion capacitance'''
        
        self.C_pn_diff_n = 0
        '''n-side small-signal diffusion capacitance per'''
        
        self.C_pn_diff_p = 0
        '''p-side small-signal diffusion capacitance per unit area'''
        
        self.C_pn_diff = 0
        '''small-signal diffusion capacitance per unit area'''
        
        self.C_pn = 0
        '''small-signal capacitance per unit area'''
        
    def dope_p_and_n(self, Na=0, Nd=0):
        self.p.Na = Na
        self.p.Nd = 0
        self.n.Na = 0
        self.n.Nd = Nd
        self.vbi = Vbi(Na, Nd)
        self.__calculate()

    def get_potential(self, side, x):
        '''gets the potential at x away from the center of the depletion region on the side specified by side.'''
        if x < -self.w_depl_n:
            return self.v_bias
        elif x > self.w_depl_p:
            return 0
        if (side == 'n'):
            return self.v_bias - v_n_depletion_edge(self.n.Nd, self.w_depl, x, self.v_bias)
        else:
            return self.v_bias - v_p_depletion_edge(self.p.Na, self.w_depl, x, self.v_bias)
        
    def bias(self, Vpn):
        self.v_bias = Vpn
        self.__calculate()
        
    def set_area(self, A):
        self.area = A
        self.__calculate()
        
    def __calculate(self):
        self.__get_depl_widths()
        self.__get_edge_potentials()
        self.__get_cur_densities()
        self.__get_e_bi_max()
        
        Nb = min(self.p.Na, self.n.Nd)
        self.V_breakdown = 60 * (cs.Eg_si/1.1242)**(3/2) * (1e16/Nb)**(3/4)
        self.Q_pn_dep = q * self.n.Nd * self.w_depl_n
        self.Q_pn_diff_n = (q*ni**2*self.n.p_diffusion_length/self.n.Nd) * (exp(q*self.v_bias/KbT) - 1)
        self.Q_pn_diff_p = -(q*ni**2*self.p.n_diffusion_length/self.p.Na) * (exp(q*self.v_bias/KbT) - 1)
        self.__get_small_signal_capacitance()
    
    def __get_depl_widths(self):
        self.w_depl = W_depl(self.p.Na, self.n.Nd, self.v_bias)
        self.w_depl_n = W_depl_n(self.p.Na, self.n.Nd, self.v_bias)
        self.w_depl_p = W_depl_p(self.p.Na, self.n.Nd, self.v_bias)
        self.depl_edge_n_conc = depl_min_conc_at_edge(self.p.n, self.v_bias)
        self.depl_edge_p_conc = depl_min_conc_at_edge(self.n.p, self.v_bias)
        
    def __get_cur_densities(self):
        self.depl_edge_n_J = electron_cur_density_at_p_depl_edge(self.p.Na, self.p.Dn, self.p.n_diffusion_length, self.v_bias)
        self.depl_edge_p_J = hole_cur_desnsity_at_n_depl_edge(self.n.Nd, self.n.Dp, self.n.p_diffusion_length, self.v_bias)
        
        self.t_rec_scr = (self.p.rec_life_n + self.n.rec_life_p)/2
        self.t_gen_scr = (self.p.gen_life_n + self.n.gen_life_p)/2
        self.Rmax_rec = ni/(2*self.t_rec_scr) * exp((q * self.v_bias)/(2*KbT) - 1)
        self.J_d_scr = q*ni*self.w_depl/(2*self.t_rec_scr) * exp((q * self.v_bias)/(2*KbT) - 1)
        self.J_d_scg = self.J_d_scr
        
        self.J_S_diff = q * self.p.Dn * ni**2 / (self.p.n_diffusion_length * self.p.Na) + q * self.n.Dp * ni**2 / (self.n.p_diffusion_length * self.n.Nd)
        self.J_S_scr = q*ni*self.w_depl / (2*self.t_rec_scr)
        
        self.J_D = self.J_S_diff * (exp(q*self.v_bias/KbT) - 1) + self.J_S_scr * (exp(q*self.v_bias/(2*KbT)) - 1)
        
        self.I_D = self.J_D * self.area
        
    def __get_edge_potentials(self):
        self.vbi_n = v_n_depletion_edge(self.n.Nd, self.w_depl_n, self.v_bias)
        self.vbi_p = v_p_depletion_edge(self.p.Na, self.w_depl_p, self.v_bias)
        
    def __get_e_bi_max(self):
        self.e_bi_max = E_p_depletion(self.p.Na, self.w_depl_p, 0)
        
    def __get_small_signal_capacitance(self):
        p1 = q*eps_si/2
        p2 = self.p.Na * self.n.Nd / (self.p.Na + self.n.Nd)
        p3 = 1/(self.vbi - self.v_bias)
        self.C_pn_dep = self.area * sqrt(p1 * p2 * p3)
        
        self.C_pn_diff_n = self.area * (q**2 * ni**2 * self.n.p_diffusion_length / (KbT * self.n.Nd)) * exp(q*self.v_bias/KbT)
        self.C_pn_diff_p = self.area * (q**2 * ni**2 * self.p.n_diffusion_length / (KbT * self.p.Na)) * exp(q*self.v_bias/KbT)
        self.C_pn_diff = self.C_pn_diff_n + self.C_pn_diff_p
    
        self.C_pn = self.C_pn_dep + self.C_pn_diff
        