from solvers.mosfet import Mosfet

class MosAmp:
    
    def __init__(self):
        self._docalc = False
        
        self.mf = Mosfet()
        self.R1 = 1
        self.R2 = 1
        self.RL = 1
        '''The load resistance of the amplifier'''
        self.RD = 1
        '''The drain resistance of the amplifier'''
        self.RS = 1
        '''The source resistance of the amplifier'''
        self.Rin = 1
        '''The input resistance of the amplifier'''
        self.Rsig = 1
        '''The signal resistance of the amplifier'''
        self.Rout = 1
        '''The output resistance of the amplifier'''
        self.C1 = 0
        self.C2 = 0
        self.C3 = 0
        self.Av = 0
        '''The voltage gain of the amplifier'''
        self.Ai = 0
        '''The current gain of the amplifier'''
        
        self.Vdd = 1
        '''The input voltage, units = V'''
        self.Vgg = 0
    
        
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
        if self._docalc and __name != '_docalc':
            self._docalc = False
            self.__calculations()
            self._docalc = True
    
    def __calculations(self):
        self.Vgg = self.Vdd * self.R1 / (self.R1 + self.R2)
        Rg = self.R1 * self.R2 / (self.R1 + self.R2)
        
        self.Rin = Rg
        
        Rp = 1/(1/self.mf.rds_sat + 1/self.RD + 1/self.RL)
        
        self.Av = -self.mf.gm * Rp * (Rg) / (self.Rsig * Rg)
        
        self.Ai = self.Av * (self.Rsig + Rg) / (self.RL)
        
        self.Rout = 1/(1/self.mf.rds_sat + 1/self.RD)
