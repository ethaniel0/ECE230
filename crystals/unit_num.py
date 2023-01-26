import math
import sympy as sp

class UnitNum:
    UNITS = {
        'm': 1,
        'ang': 1e-10,
        'Å': 1e-10,
        'angstroms': 1e-10,
        'cm': 1e-2,
        'mm': 1e-3,
        'nm': 1e-9,
        'pm': 1e-12,
        'fm': 1e-15,
        'rad': 1,
        'deg': sp.pi / 180,
    }
    BASES = {
        'm': 'm',
        'ang': 'm',
        'Å': 'm',
        'angstroms': 'm',
        'cm': 'm',
        'mm': 'm',
        'nm': 'm',
        'pm': 'm',
        'fm': 'm',
        'rad': 'rad',
        'deg': 'rad'
    }
    def __init__(self, num, units: str, convert=False) -> None:
        if isinstance(num, UnitNum):
            if convert:
                self = num.to(units)
            else:
                self.num = num.num
                self.units = units
        else:
            self.num = num
            self.units = units

    def to(self, unit: str):
        if unit == self.units or unit not in UnitNum.UNITS:
            return self.copy()
        if UnitNum.BASES[unit] != UnitNum.BASES[self.units]:
            return self.copy()
        
        num = self.num * UnitNum.UNITS[self.units] / UnitNum.UNITS[unit]
        return UnitNum(num, unit)
    
    def copy(self):
        return UnitNum(self.num, self.units)

    def __str__(self):
        return f'{self.num} {self.units}'
    def __repr__(self):
        return f'{self.num} {self.units}'
    def __add__(self, other):
        if self.units == other.units or not isinstance(other, UnitNum):
            return UnitNum(self.num + other.num, self.units)
        elif isinstance(other, UnitNum) and UnitNum.BASES[self.units] == UnitNum.BASES[other.units]:
            return UnitNum(self.num + other.to(self.units).num, self.units)
        else:
            raise ValueError(f'Cannot add units of {self.units} and {other.units}')
    def __sub__(self, other):
        if self.units == other.units:
            return UnitNum(self.num - other.num, self.units)
        elif isinstance(other, UnitNum) and UnitNum.BASES[self.units] == UnitNum.BASES[other.units]:
            return UnitNum(self.num - other.to(self.units).num, self.units)
        else:
            raise ValueError(f'Cannot subtract units of {self.units} and {other.units}')
    def __mul__(self, other):
        if isinstance(other, UnitNum):
            return UnitNum(self.num * other.num, f'{self.units}*{other.units}')
        else:
            return UnitNum(self.num * other, self.units)
    def __truediv__(self, other):
        if isinstance(other, UnitNum):
            if self.units == other.units:
                return self.num / other.num
            return UnitNum(self.num / other.num, f'{self.units} / {other.units}')
        else:
            return UnitNum(self.num / other, self.units)
    def __pow__(self, other):
        if isinstance(other, UnitNum):
            return UnitNum(self.num ** other.num, f'{self.units} ^ {other.units}')
        else:
            return UnitNum(self.num ** other, self.units)
    def __eq__(self, other):
        if isinstance(other, UnitNum):
            return self.num == other.num and self.units == other.units
        else:
            return self.num == other
    def __ne__(self, other):
        if isinstance(other, UnitNum):
            return self.num != other.num or self.units != other.units
        else:
            return self.num != other
    def __lt__(self, other):
        if isinstance(other, UnitNum):
            if self.units == other.units:
                return self.num < other.num
            else:
                raise ValueError(f'Cannot compare units of {self.units} and {other.units}')
        else:
            return self.num < other
    def __le__(self, other):
        if isinstance(other, UnitNum):
            if self.units == other.units:
                return self.num <= other.num
            else:
                raise ValueError(f'Cannot compare units of {self.units} and {other.units}')
        else:
            return self.num <= other
    def __gt__(self, other):
        if isinstance(other, UnitNum):
            if self.units == other.units:
                return self.num > other.num
            else:
                raise ValueError(f'Cannot compare units of {self.units} and {other.units}')
        else:
            return self.num > other
    def __ge__(self, other):
        if isinstance(other, UnitNum):
            if self.units == other.units:
                return self.num >= other.num
            else:
                raise ValueError(f'Cannot compare units of {self.units} and {other.units}')
        else:
            return self.num >= other
    def __abs__(self):
        return UnitNum(abs(self.num), self.units)
    def __neg__(self):
        return UnitNum(-self.num, self.units)
    def __pos__(self):
        return UnitNum(+self.num, self.units)
    def __round__(self, ndigits=None):
        return UnitNum(round(self.num, ndigits), self.units)
    def __floor__(self):
        return UnitNum(math.floor(self.num), self.units)
    def __ceil__(self):
        return UnitNum(math.ceil(self.num), self.units)
    def __trunc__(self):
        return UnitNum(math.trunc(self.num), self.units)
    def __float__(self):
        return float(self.num)
    def __int__(self):
        return int(self.num)
    def __hash__(self):
        return hash(self.num)
    def __format__(self, format_spec):
        return format(self.num, format_spec)
    def __radd__(self, other):
        return self + other
    def __rsub__(self, other):
        return -self + other
    def __rmul__(self, other):
        return self * other
    def __rtruediv__(self, other):
        return other / self.num
    def __rpow__(self, other):
        return other ** self.num
    def __iadd__(self, other):
        self = self + other
        return self
    def __isub__(self, other):
        self = self - other
        return self
    def __imul__(self, other):
        self = self * other
        return self
    def __itruediv__(self, other):
        self = self / other
        return self
    def __ipow__(self, other):
        self = self ** other
        return self
    def __matmul__(self, other):
        return self * other
    def __rmatmul__(self, other):
        return self * other
    def __imatmul__(self, other):
        self = self @ other
        return self

if __name__ == '__main__':
    u1 = UnitNum(1, 'm')
    u2 = UnitNum(2, 'm')
    print(u1 / u2)