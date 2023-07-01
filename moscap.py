from sympy import symbols, solve, ln, sqrt, exp
from solvers.constants import KbToq, ni, eps_si, q, KbT, eps_ox
from solvers.silicon import Silicon

def Vgb(Vox=symbols("Vox"), Vscb=symbols("Vscb"), phi_pm=symbols("phi_pm")):
    '''calculates the gate-to-bulk bias voltage'''
    return Vox + Vscb + phi_pm

def charge_density_conservation(solveFor, Qg=symbols("Q_g"), Qscb=symbols("Q_scb"), Qit=symbols("Q_it"), Qf=symbols("Q_f")):
    '''calculates the charge density in the substrate region
    solveFor = the variable to solve for (Qscb, Qit, Qf, Qg)
    '''
    solveFor = symbols(solveFor)
    return solve(Qg + Qscb + Qit + Qf, solveFor)[0]

def gate_charge_density(eox=symbols("ϵ_ox"), Exox=symbols("E_xox")):
    '''calculates Q_G = gate charge density'''
    return eox * Exox

def V_contact(Na_ionized):
    '''calculates phi_pm = the contact potential difference between the neutral bulk region and the aluminum back-side contact'''
    return -0.51165 - KbToq * ln(Na_ionized / ni)

def fermi_potential(Na_ionized):
    '''calculates phi_fb = the fermi potential in the substrate region'''
    return KbToq * ln(Na_ionized / ni)

def width_at_Vtn(Nab, Vscb):
    frac = 2*eps_si/(q * Nab)
    return sqrt(frac * Vscb)

def substrate_debeye_length(Nab):
    '''calculates the substrate debeye length'''
    return sqrt(eps_si * KbT / (q**2 * Nab))

def calc_space_charge_density(Vscb, Ldb, nb, pb):
    '''calculates Q_SC,B = the space-charge density in the p-type silicon substrate at the onset of strong inversion'''
    beta = q/KbT
    p1 = sqrt(2)*eps_si * KbToq / Ldb
    inner_1 = exp(-beta * Vscb) + beta*Vscb - 1
    inner_2 =  exp(beta * Vscb) - beta*Vscb - 1
    return -p1 * sqrt(inner_1 + (nb/pb) * inner_2)

def calc_e_xox(Qg):
    return Qg / symbols("ϵ_ox")

class AluminumMoscap(Silicon):
    '''
    An aluminum MOSCAP.
    '''
    def __init__(self):
        super().__init__()
        
        self._docalc = False
        
        # MOSCAP DESIGN PARAMETERS
        self.Nab = symbols("N_ab")
        '''the acceptor concentration in the substrate region. units = cm^-3'''
        self.Xox = symbols("X_ox")
        '''the gate-oxide thickness. units = cm'''
        self.Xw = symbols("X_w")
        '''the wafer thickness. units = cm'''
        self.Qf = symbols("Q_f")
        '''the fixed-oxide charge density.'''
        self.phi_pm = symbols("phi_pm")
        '''the contact-potential difference between the p-type substrate and the backside metal contact. units = V'''
        self.phi_fb = symbols("phi_fb")
        '''the fermi potential in the substrate region. units = V'''
        
        self.VGB = symbols("V_GB")
        '''V_GB = the gate-to-bulk bias voltage (the applied voltage). units = V'''
        
        # BIAS-DEPENDENT VARIABLES
        self.Qg = symbols("Q_g")
        '''the charge density in the metal-gate region. units = C/cm^2'''
        self.Qit = symbols("Q_it")
        '''the interface-trap charge density. units = C/cm^2'''
        self.Qscb = symbols("Q_scb")
        '''the charge density in the substrate space-charge region. units = C/cm^2'''
        self.QinvB = symbols("Q_invB")
        '''the substrate inversion charge density. units = C/cm^2'''
        
        self.E_gate_oxide = symbols("E_xox")
        '''E_x,OX = the transverse electric field in the gate-oxide layer if it is assumed to be charge-free. units = V/cm'''
        self.V_gate_oxide = symbols("Vox")
        '''V_OX = the voltage across the gate-oxide layer. units = V'''
        self.Vscb = symbols("Vscb")
        '''V_SC,B = the voltage across the substrate space-charge region. units = V'''
        self.v_space_charge_at_Vtn = symbols("v_SCB(V_TN)")
        '''v_SC,B(V_TN) = the threshold voltage across the substrate space-charge region. units = V'''
        
        # BIAS-DEPENDENT MOS-VARIABLE DISTRIBUTIONS
        self.gate_oxide_potential = symbols("phi_ox")
        '''phi_ox, electrostatic-potential distribution in the gate-oxide layer. units = V'''
        self.E_space_charge = symbols("E_x_B")
        '''E_x,B the transverse electric-field distribution in the substrate space-charge region'''
        self.n_space_charge = symbols("n_B")
        '''n_B = electron concentration distribution in the substrate space-charge region'''
        self.p_space_charge = symbols("p_B")
        '''p_B = hole concentration distribution in the substrate space-charge region'''
        self.net_charge_dist_space_charge = symbols("varphi_B")
        '''varphi_B = net charge distribution in the substrate space-charge region'''
        
        self.debye_length = symbols("L_BD")
        '''L_B,D = the substrate debeye length'''
        
        self.Vtn = symbols("V_TN")
        '''V_TN = the threshold voltage'''
        
        self.VFB = symbols("V_FB")
        '''V_FB = the flat-band voltage'''
        
        self.Cox = symbols("C_ox")
        '''C_ox = the gate-oxide capacitance'''
        
        self.Cgb = symbols("C_gb")
        '''C_gb = the gate-substrate capacitance. In the inversion range, this is the low-frequency value.'''
        
        self.Cgb_HF = symbols("C_gb")
        '''C_gb = the gate-substrate capacitance'''
        
        self.Cscb = 0
        '''C_SC,B = the substrate space-charge region capacitance. In the inversion range, this is the low-frequency value.'''
        
        self.Cscb_HF = 0
        '''C_SC,B-HF = the substrate space-charge region capacitance for high frequencies'''
        
        self.Cgb_threshold = symbols("C_gb")
        '''C_gb = the gate-substrate capacitance'''
        
        self.w_dep = symbols("w_dep")
        '''Width of the depletion region'''
        
        self.last_modified = None
        self.value_lock = []
        
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
        super().__setattr__(__name, __value)
        self.__dict__[__name] = __value
        if __name == 'Nab':
            self.Na = self.Nab
            self.Nd = 0
        if self._docalc and __name != '_docalc':
            self._docalc = False
            self.__calculations()
            self._docalc = True
        
    def __calculations(self):
        self.phi_pm = V_contact(self.Nab)
        self.phi_fb = fermi_potential(self.Nab)
        self.Cox = eps_ox/self.Xox
        self.VFB = self.phi_pm - (self.Qf + self.Qit)/self.Cox
        self.Vtn = self.VFB + 2*self.phi_fb + sqrt(4*q*eps_si*self.Nab*self.phi_fb)/self.Cox
        
        self.v_space_charge_at_Vtn = 2*self.phi_fb
        self.width_at_Vtn = width_at_Vtn(self.Nab, self.v_space_charge_at_Vtn)
        self.debye_length = substrate_debeye_length(self.Nab)
        
        # default accumulation range if VGB isn't defined
        if not self.isnumber(self.VGB) or self.VGB < self.VFB:
            self.V_gate_oxide = self.VGB - self.phi_pm
            self.E_gate_oxide = self.V_gate_oxide/self.Xox
            self.Qg = eps_ox * self.E_gate_oxide
            self.Qscb = -self.Cox
            self.Vscb = self.VGB - self.phi_pm + (self.Qf + self.Qit)/self.Cox + self.Qscb/self.Cox
            self.w_dep = 0
            self.Cgb = self.Cox * self.area
        
        # at the flat band voltage
        elif self.VGB == self.VFB:
            self.Qscb = 0
            self.Vscb = 0
            self.Qg = -self.Qf - self.Qit
            self.E_gate_oxide = self.Qg/eps_ox
            self.V_gate_oxide = self.Xox * self.E_gate_oxide
            self.w_dep = 0
            self.Cgb = self.Cox * self.area
            
        # depletion range
        elif self.VFB < self.VGB < self.Vtn:
            # print('in depl range', self.VGB)
            self.Vscb = self.VGB - self.VFB \
                        + q*eps_si*self.Nab/(self.Cox**2) * \
                        (1 - sqrt(1 + 2*self.Cox**2/(q*eps_si*self.Nab) * (self.VGB - self.VFB)))
            self.E_gate_oxide = sqrt(2*q*eps_si*self.Nab*self.Vscb)/eps_ox - (self.Qf + self.Qit)/eps_ox
            self.V_gate_oxide = self.Xox * self.E_gate_oxide
            self.Qg = sqrt(2*q*eps_si*self.Nab*self.Vscb) - (self.Qf + self.Qit)
            self.Qscb = -(self.Qg + self.Qf + self.Qit)
            self.w_dep = sqrt(2*eps_si/(q*self.Nab) * self.Vscb)
            
            self.Cscb = eps_si/self.w_dep
            self.Cgb = self.Cox * self.Cscb / (self.Cox + self.Cscb)
            self.Cscb *= self.area
            self.Cgb *= self.area
            
        
        else: # inversion range
            self.Vscb = 2*self.phi_fb
            self.w_dep = self.width_at_Vtn
            self.Vox = self.VGB - 2*self.phi_fb - self.phi_pm
            self.E_gate_oxide = self.Vox/self.Xox
            self.V_gate_oxide = self.Xox * self.E_gate_oxide
            self.Qg = self.E_gate_oxide * eps_ox
            self.Qscb = -(self.Qg + self.Qf + self.Qit)
            self.QinvB = -self.Cox * (self.VGB - self.Vtn)
            
            self.Cscb_HF = eps_si/self.w_dep
            self.Cgb_HF = self.Cox * self.Cscb / (self.Cox + self.Cscb)
            self.Cscb_HF *= self.area
            self.Cgb_HF *= self.area
            
            self.Cgb = self.Cox * self.area
            
        self.E_space_charge = -self.Qscb/eps_si
            
        
    
    def isnumber(self, value):
        try:
            if value > 12:
                pass
            return True
        except:
            return False
        