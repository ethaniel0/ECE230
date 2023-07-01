from sympy import symbols, sqrt, Rational, ln
import sympy as sp
from math import pi
from solvers import moscap
from solvers import constants as cs

class Mosfet(moscap.AluminumMoscap):
    def __init__(self):
        super().__init__()
        
        self._docalc = False
        
        self.channel_length = symbols("L_ch")
        '''L_ch, the channel length of the MOSFET. units = cm'''
        self.channel_width = symbols("W_ch")
        '''W_ch, the channel width of the MOSFET units = cm'''
        self.trap_charge_densities = symbols("Qf+Qit")
        
        '''Qf+Qit, the fixed-oxide and interface-trap charge densities'''
        self.Vgs = symbols("V_GS")
        '''V_GS, the gate-source voltage, units = V'''
        self.Vgd = symbols("V_GD")
        '''V_GD, the gate-drain voltage, units = V'''
        self.Vds = 0
        '''V_DS, the drain-to-source voltage, units = V'''
        self.Vdsat = symbols("V_DSat")
        '''V_DSat, the drain-to-source saturation voltage, units = V'''
        self.Vbs = symbols("V_BS")
        '''V_BS, the bulk-to-source voltage, units = V'''
        self.Vfb_n0 = symbols("V_fb_n0")
        '''V_fb_n0, the flatband voltage at zero gate voltage'''
        
        self.mu_n_ch = symbols("mu_n_ch")
        '''mu_n_ch, the electron mobility in the channel region'''
        
        self.Id = symbols("I_Dsat")
        '''I_D, the drain current'''
        
        self.ch_len_modulation = 0
        '''λ_N channel length modulation factor'''
       
        
        self.gm = symbols("g_m")
        '''g_m, the forward transconductance, units = S'''
        self.gds = symbols("g_ds")
        '''g_ds, the drain-source conductance, units = S'''
        self.Q_dep_B = symbols("Q_dep_B")
        '''Q_dep_B, the substrate depletion charge density'''
        self.Q_inv_B = symbols("Q_inv_B")
        '''Q_inv_B, the substrate inversion charge density'''
        
        self.rds_sat = symbols("r_ds_sat")
        '''r_ds_sat, the saturation drain-source resistance'''
        
        self.Kn = symbols("K_n")
        '''K_n = µ_n,ch*C_ox*W_ch/2L_ch, the electron affinity of the substrate. units = A/V^2'''
        
        self.K_N = symbols("K_N")
        '''K_N = µ_n,ch*C_ox*W_ch/L_ch, the electron affinity of the substrate. units = A/V^2'''
        
        self.bigQg = symbols("Q_G")
        '''element of gate charge between y=0 and y=L_ch'''
        
        self.bigQg_lin = symbols("Q_G_lin")
        '''element of gate charge between y=0 and y=L_ch in the linear range'''
        
        self.total_Cox = symbols("C_ox_tot")
        '''C_ox_tot, the total oxide capacitance = Cox * W_ch * L_ch'''
        
        self.C_gs = symbols("C_gs")
        '''C_gs, the gate-to-source capacitance'''
        self.C_gs_overlap = symbols("C_gs)ov")
        '''C_gs, the gate-to-source overlap capacitance'''
        self.C_gd = symbols("C_gd")
        '''C_gd, the gate-to-drain capacitance'''
        self.C_gd_overlap = symbols("C_gd_ov")
        '''C_gd_ov, the gate-to-drain overlap capacitance'''
        
        self.Lgsov = symbols("L_gs_ov")
        '''L_gs_ov, the gate-to-source overlap length'''
        self.Lgdov = symbols("L_gd_ov")
        '''L_gd_ov, the gate-to-drain overlap length'''
        
        self.bulk_body_effect_coeff = symbols("gamma")
        '''γ, the bulk-body effect coefficient'''
        
        self.vbis = []
        '''V_bi for each interface: channel/drain, side-wall/drain, side-wall/drain, side-wall/drain, substrate/drain'''
        
        self.ft = symbols("f_t")
        
        self.last_modified = None
        self.value_lock = []
        
        self._docalc = True
        
    def __calculate(self):
        self._docalc = False
        
        self.Q_dep_B = self.depl_charge_density(self.Nab, self.phi_fb, self.Vbs)
        self.Kn = self.mu_n_ch * self.Cox * self.channel_width / (2 * self.channel_length)
        self.K_N = 2*self.Kn
        
        if (self.isnumber(self.K_N) and not self.isnumber(self.Kn)):
            self.Kn = self.K_N / 2
        
        self.Vtn = Mosfet.thresh_voltage(self.VFB, self.phi_fb, self.Cox, self.Nab, self.Vbs)
        self.bulk_body_effect_coeff = sqrt(2*cs.q*cs.eps_si*self.Nab)/self.Cox
        self.Q_inv_B = self.inv_charge_density_easy(self.Cox, self.Vgs, self.Vtn)
        self.Vdsat = self.Vgs - self.Vtn
        self.Vgd = self.Vgs - self.Vds
        self.total_Cox = self.Cox * self.channel_width * self.channel_length
        self.bigQg_lin = self.lin_gate_charge(self.Cox * self.channel_width * self.channel_length, self.Vgs, self.Vds, self.Vtn)
        
        self.__set_vbis()
        
        self.ft = (3*self.mu_n_ch) / (4*pi*self.channel_length**2) * (self.Vgs - self.Vtn)
        
        vgs_sym = symbols("Vgs")
        
        if self.isnumber(self.Vgs) and self.isnumber(self.Vtn) and self.Vgs < self.Vtn:
            current = 0
        elif self.isnumber(self.Vds) and self.isnumber(self.Vdsat) and self.Vds >= self.Vdsat:
            current = self.drain_current_sat_kn(self.Kn, vgs_sym, self.Vtn, self.ch_len_modulation, self.Vds, self.Vdsat)
            self.C_gs = (2/3) * self.total_Cox
            self.C_gd = 0
        else:
            current = self.drain_current_lin_kn(self.Kn, vgs_sym, self.Vtn, self.Vds)
            self.C_gs = (2/3) * self.total_Cox * (1 - (self.Vgd - self.Vtn)**2 / (self.Vgs + self.Vgd - 2*self.Vtn)**2 )
            self.C_gd = (2/3) * self.total_Cox * (1 - (self.Vgs - self.Vtn)**2 / (self.Vgs + self.Vgd - 2*self.Vtn)**2 )
        
        self.C_gs_overlap = cs.eps_ox*self.channel_width*self.Lgsov/self.Xox
        self.C_gd_overlap = cs.eps_ox*self.channel_width*self.Lgdov/self.Xox
            
        self.Id = 0 if current == 0 else current.subs(vgs_sym, self.Vgs)
        
        id_diff = 0 if current == 0 else current.diff(vgs_sym)
        self.gm = 0 if current == 0 else id_diff.subs(vgs_sym, self.Vgs)
        self.gds = 0 if current == 0 else self.mu_n_ch * self.Cox * self.channel_width/self.channel_length * (self.Vgs - self.Vtn)**2 * self.ch_len_modulation
        self.rds_sat = sp.oo if self.gds == 0 else 1/self.gds
        
        Idlin = self.drain_current_lin_kn(self.Kn, self.Vgs, self.Vtn, self.Vds)
        self.bigQg = self.gate_charge(self.mu_n_ch, self.channel_width, self.Cox, Idlin, self.Vgs, self.Vtn, self.Vds)
        
        self._docalc = True
        
    
    def lock(self, name):
        self.value_lock.append(name)
        
    def set_lock(self, name, value):
        self.__setattr__(name, value)
        self.lock(name)
        
    def unlock(self, name):
        if name in self.value_lock:
            self.value_lock.remove(name)
        
    def __setattr__(self, __name: str, __value) -> None:
        if 'value_lock' in self.__dict__ and __name in self.value_lock:
            return
        
        if '_docalc' in self.__dict__ and 'last_modified' in self.__dict__ and not self._docalc and __name == self.last_modified:
            return
        
        super().__setattr__(__name, __value)
        self.__dict__[__name] = __value
        
        if __name == 'trap_charge_densities':
            self.Qf = self.trap_charge_densities
            self.Qit = 0
        
        if self._docalc and __name != '_docalc':
            self._docalc = False
            self.last_modified = __name
            self.__calculate()
            self._docalc = True
    
    def __set_vbis(self):
        vbi1 = cs.KbToq * ln(self.Nab * self.Nd/cs.ni**2)
        vbi2 = cs.KbToq * ln(self.Na * self.Na/cs.ni**2)
        self.vbis = [vbi1, vbi2, vbi2, vbi2, vbi1]
            
    @staticmethod
    def thresh_voltage(Vfb_n0, phi_fb, Cox, Nab, Vbs):
        return Vfb_n0 + 2*phi_fb + (1/Cox)*sqrt( 2 * cs.q * cs.eps_si * Nab * (2*phi_fb - Vbs) )
        
    @staticmethod
    def depl_charge_density(Nab, phi_fb, Vbs):
        '''
        MOSFET Level-1 Model for Q_DEP,B
        '''
        return -sqrt( 2 * cs.q * cs.eps_si * Nab * (2*phi_fb - Vbs) )
    
    @staticmethod
    def inv_charge_density(Cox, Vgs, Vfbn0, phi_fb, Nab, Vbs):
        '''
        MOSFET Level-1 Model for Q_INV,B
        '''
        return sqrt( 2 * cs.q * cs.eps_si * Nab * (2*phi_fb - Vbs) ) - (Vgs - Vfbn0 - 2*phi_fb) / Cox
    
    @staticmethod
    def inv_charge_density_easy(Cox, Vgs, Vtn):
        return -Cox * (Vgs - Vtn)

    @staticmethod
    def drain_current_sat(mu_n_ch, Cox, Wch, Lch, Vgs, Vtn, lambda_n=0, Vds=0, Vdsat=0):
        return mu_n_ch * Cox * (Wch / (2 * Lch)) * (Vgs - Vtn)**2 * (1 + lambda_n * (Vds - Vdsat))
    
    @staticmethod
    def drain_current_sat_kn(Kn, Vgs, Vtn, lambda_n=0, Vds=0, Vdsat=0):
        return Kn * (Vgs - Vtn)**2 * (1 + lambda_n * (Vds - Vdsat))
    
    @staticmethod
    def drain_current_lin(Kn, Vgs, Vtn, Vds):
        return 2*Kn * ( (Vgs - Vtn)*Vds - 0.5*Vds**2 )
    
    @staticmethod
    def drain_current_lin_kn(Kn, Vgs, Vtn, Vds):
        return 2 * Kn * ((Vgs - Vtn)*Vds - 0.5*Vds**2 )
    
    @staticmethod
    def gate_charge(m_unch, Wch, Cox, Idlin, Vgs, Vtn, Vds):
        return (m_unch * Wch**2 * Cox**2 / Idlin) * ( (Vgs-Vtn)**2*Vds - (Vgs-Vtn)*Vds**2 + (1/3)*Vds**3 )
    
    @staticmethod
    def lin_gate_charge(total_Cox, Vgs, Vds, Vtn):
        Vgd = Vgs - Vds
        if (Vgd - Vtn)**2 == (Vgs - Vtn)**2:
            return 0
        return (2/3) * total_Cox * ( (Vgd - Vtn)**3 - (Vgs - Vtn)**3 ) / ( (Vgd - Vtn)**2 - (Vgs - Vtn)**2 )
    
    def isnumber(self, value):
        try:
            if value > 12:
                pass
            return True
        except:
            return False
