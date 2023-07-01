from solvers.mosfet import Mosfet

class MosInverter:
    
    def __init__(self):
        self._docalc = False
        
        self.mf = Mosfet()
        self.RL = 1
        self.Iout = 0
        self.ID = 0
        self.VDD = 0
        
        self.Vout = 0
        self.Vin = 0
    
        
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
        p1 = self.VDD/self.RL - self.mf.K_N/2 * (self.Vin - self.mf.Vtn) * (1 - self.mf.ch_len_modulation * (self.Vin - self.mf.Vtn))
        p2 = self.mf.ch_len_modulation * self.mf.K_N/2 * (self.Vin - self.mf.Vtn)**2 + 1/self.RL
        self.Vout = p1/p2
        
        self.Id = (self.VDD - self.Vout)/self.RL
        