from solvers import carriers as eqs
from solvers import constants as cs
from sympy import sqrt, symbols

class Silicon:
    def __init__(self):
        self._docalc = False
        
        self.ATOM_CONCENTRATION = 5e22
        '''5e22 atoms/cm^3' in a silicon crystal'''
        self.UNIT_CELL_VOLUME = 5.43095**3
        '''Volume of the unit cell of silicon = 5.43095 Å^3'''
        self.LATTICE_CONSTANT = 5.43095
        '''Lattice constant of silicon = side length of Si unit cell =  5.43095 Å'''
        self.ATOMIC_MASS=28.0855
        '''Atomic mass of Silicon = 28.0855 g/mol'''
        self.BANDGAP = 1.1242
        '''Bandgap of Silicon = 1.1242 eV'''

        self.is_n_type: bool = False
        self.is_p_type: bool = False
        
        self.Na = symbols("Na") 
        '''acceptor doping concentration. units = cm^-3'''
        
        self.Nd = symbols("Nd") 
        '''donor doping concentration. units = cm^-3'''
        
        self.compensation_threshold = 5
        
        self.n = symbols("n") 
        '''electron concentration. units = cm^-3'''
        self.p = symbols("p") 
        '''hole concentration. units = cm^-3'''
        
        self.mu_p = symbols("µ_p") 
        '''μ_p, hole mobility. units = cm^2/V⋅s'''
        self.mu_n = symbols("µ_n") 
        '''μ_n, electron mobility units = cm^2/V⋅s'''
        
        self.e_field = 0
        '''electric field strength units = V/cm'''
        
        self.rec_life_n = symbols("τ_rec_N") 
        '''τ_rec,P = τ_n,P - electron recombination lifetime. units = s'''
        self.rec_life_p = symbols("τ_rec_P") 
        '''τ_rec,N = τ_p,N - hole recombination lifetime. units = s'''
        
        self.rec_len_n = symbols("L_rec_N")
        '''L_rec,P = L_n,P - electron recombination length. units = cm'''
        self.rec_life_p = symbols("L_rec_P")
        '''L_rec,N = L_p,N - hole recombination length. units = cm'''
        
        self.gen_life_n = symbols("τ_gen_N") 
        '''τ_gen,P = τ_n,P - minority electron generation lifetime. units = s'''
        self.gen_life_p = symbols("τ_gen_P") 
        '''gen,N = τ_p,N - minority hole generation lifetime. units = s'''
        
        self.gen_len_n = symbols("L_gen_N")
        self.gen_len_p = symbols("L_gen_P")
        
        self.n_diffusion_length = symbols("L_diff_n")
        '''L_n,P = L_n - electron diffusion length. units = cm'''
        self.p_diffusion_length = symbols("L_diff_p") 
        '''L_p,N = L_p - hole diffusion length. units = cm'''
        
        self.n_debye_length = symbols("L_DN") 
        ''''L_D,N = debye length in n-type region. units = cm'''
        self.p_debye_length = symbols("L_DP") 
        '''L_D,P = debye length in p-type region. units = cm'''
        
        self.Dn = symbols("Dn") 
        '''D_n - electron diffusivity. units = cm^2/s'''
        self.Dp = symbols("Dp")
        '''D_p - hole diffusivity. units = cm^2/s'''
        
        self.conductivity = symbols("σ") 
        '''σ - electrical conductivity. units = S/cm'''
        self.resistivity = symbols("ρ") 
        '''ρ - electrical resistivity. units = Ω⋅cm'''
        self.resistance = symbols("R") 
        '''R - electrical resistance. units = Ω'''
        self.dielectric_relxation_time = 1
        '''τ_diel - dielectric relaxation time. units = s'''
        
        self.length = 1
        '''length of the sample. units = cm'''
        self.area = 1
        '''area of the sample. units = cm^2'''
        
        self.v_n_drift = symbols("v_n_drift")
        '''electron drift velocity. units = cm/s'''
        self.J_n_drift = symbols("J_n_drift") 
        '''electron drift current density. units = A/cm^2'''
        
        self.v_p_drift = symbols("v_p_drift")
        '''hole drift velocity. units = cm/s'''
        self.J_p_drift = symbols("J_p_drift")
        '''hole drift current density. units = A/cm^2'''
        
        self.J_x_drift = None
        '''total drift current density. units = A/cm^2'''
        
        self.is_degenerate = False
        '''If the material is degenerate (n > N_c in n-type or p > N_v in p-type)'''
        
        self.last_modified = None
        self.value_lock = []
        
        self._docalc = True
        
        self.__calculate()
    
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
        self.__dict__[__name] = __value
        if self._docalc and __name != '_docalc':
            self._docalc = False
            self.__calculate()
            self._docalc = True
    
    def __calculate(self):
        # print('calculating from silicon')
        self._docalc = False
        
        if  (type(self.Na) in [int, float] and type(self.Nd) in [int, float]):
        
            self.n = self.Nd
            self.p = self.Na
                
            if self.Na > self.Nd:
                self.is_p_type = True
                self.is_n_type = False
            elif self.Nd > self.Na:
                self.is_n_type = True
                self.is_p_type = False
            
        self.__calc_n_p()
        self.__calc_mobilty()
        self.__calc_rec_gen()
        self.__calc_diffusivity()
        self.__calc_diff_length()
        self.__calc_conductivity_resistivity()
        self.__calc_debye_length()
        self.__calc_v_j_drift()
        self.__calc_diff_length()
        
        self._docalc = True
        
    def __calc_v_j_drift(self):
        self.v_n_drift = eqs.v_drift_n(self.mu_n, self.e_field)
        self.v_p_drift = eqs.v_drift_p(self.mu_p, self.e_field)
        
        self.J_n_drift = eqs.J_drift_n(self.mu_n, self.n, self.e_field)
        self.J_p_drift = eqs.J_drift_p(self.mu_p, self.p, self.e_field)
        
        if self.conductivity and self.e_field:
            self.J_x_drift = self.conductivity*self.e_field
        else:
            self.J_x_drift = self.J_n_drift + self.J_p_drift
    
    def __calc_n_p(self):
        if not (type(self.Na) in [int, float] and type(self.Nd) in [int, float]):
            return
        # if compensated silicon
        if self.Na > 0 and self.Nd > 0:
            rel = self.Na/self.Nd
            if rel > (1/self.compensation_threshold) and rel < self.compensation_threshold:
                if self.is_p_type:
                    self.p = (1/2) * ((self.Na - self.Nd) + sqrt((self.Na - self.Nd)**2 + 4*cs.ni**2))
                    self.n = cs.ni**2/self.p
                elif self.is_n_type:
                    self.n = (1/2) * ((self.Nd - self.Na) + sqrt((self.Nd - self.Na)**2 + 4*cs.ni**2))
                    self.p = cs.ni**2/self.n
        elif self.Na > self.Nd:
            self.p = self.Na
            self.n = cs.ni**2/self.Na
        elif self.Na < self.Nd:
            self.n = self.Nd
            self.p = cs.ni**2/self.Nd
        else: # Nd == Na
            self.p = cs.ni
            self.n = cs.ni
        
        self.is_degenerate = (self.is_n_type and self.n > cs.Nc) or (self.is_p_type and self.p > cs.Nv)
    
    def __calc_rec_gen(self):
        
        if not (self.is_n_type or self.is_p_type):
            return
        
        # reset the numbers
        self.rec_life_n = None
        self.rec_life_p = None
        self.gen_life_n = None
        self.gen_life_p = None
        
        if self.is_p_type: # p-type
            self.rec_life_n = eqs.rec_life_p(self.Na)
            self.gen_life_n = eqs.gen_life_p(self.Na)
            
        else: # n-type
            self.rec_life_p = eqs.rec_life_n(self.Nd)
            self.gen_life_p = eqs.gen_life_n(self.Nd)
        
    def __calc_diff_length(self):
        '''
        this is wrong, look at equation
        '''
        if not (self.is_n_type or self.is_p_type):
            return
        
        self.p_diffusion_length = None
        self.n_diffusion_length = None
        if self.is_n_type:
            self.p_diffusion_length = eqs.diffusion_length(self.Dp, self.rec_life_p)
        elif self.is_p_type:
            self.n_diffusion_length = eqs.diffusion_length(self.Dn, self.rec_life_n)
        
    def __calc_mobilty(self):
        if not (self.is_n_type or self.is_p_type):
            return
        
        if self.is_n_type:
            self.mu_n = eqs.mu_n_maj(self.Nd)
            self.mu_p = eqs.mu_p_min(self.Nd)
        else:
            self.mu_n = eqs.mu_n_min(self.Na)
            self.mu_p = eqs.mu_p_maj(self.Na)
    
    def __calc_diffusivity(self):
        self.Dn = eqs.diffusivity(self.mu_n)
        self.Dp = eqs.diffusivity(self.mu_p)
        
        if not (self.is_n_type or self.is_p_type):
            return
        
        if self.is_n_type:
            self.rec_len_p = sqrt(self.Dp * self.rec_life_p)
            self.gen_len_p = sqrt(self.Dp * self.gen_life_p)
        if self.is_p_type:
            self.rec_len_n = sqrt(self.Dn * self.rec_life_n)
            self.gen_len_n = sqrt(self.Dn * self.gen_life_n)
    
    def __calc_conductivity_resistivity(self):
        self.conductivity = eqs.conductivity(self.mu_n, self.mu_p)
        self.resistivity = 1/self.conductivity
        self.resistance = eqs.resistance_from_resistivity(self.resistivity, self.length, self.area)
        self.dielectric_relxation_time = cs.eps_si / self.conductivity
        self.__calc_debye_length()
    
    def __calc_debye_length(self):
        if not (self.is_n_type or self.is_p_type):
            return
        
        self.n_debye_length = None
        self.p_debye_length = None
        
        if self.is_n_type:
            self.n_debye_length = eqs.debye_length(self.Dn, self.dielectric_relxation_time)
        elif self.is_p_type:
            self.p_debye_length = eqs.debye_length(self.Dp, self.dielectric_relxation_time)
    