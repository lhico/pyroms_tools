import numpy as np


class dQdT(object):
    def __init__(self, scalarwind, pressure, T, rhoair=1):
        self.T = T
        self.Ua = scalarwind
        self.rhoair = rhoair
        self.Pa = pressure


    def dQdT_infrared(self):
        """Equation 5a"""
        sigma = 5.67 * 10**-8  # W m**-2*K**-4
        self.dQdT_ir = -4 * sigma * self.T**3


    def dQdT_sensibleheat(self):
        """Equation 5b"""
        Cp = 1.0048 * 10**3  # air specific heat [J kg**-1 **K-1 ]
        Ch = 10**-3  # Bulk transfer coeficient for sensible heat
        self.dQdT_sh = self.rhoair * Cp * Ch * self.Ua


    def dQdT_latentheat(self):
        """Equation 5c"""

        Ce = 1.15e-3
        L = 2.508 * 1e6  # [J kg**-1]
        es = lambda T: 10 ** (9.4051 - 2353/T)  # saturated water pressure vapo
        qs = lambda p,T: 0.622/p * es(T)  # humidity

        self.dQdT_lh = -self.rhoair * Ce * L * self.Ua * \
                       2353*np.log(10)*qs(self.Pa, self.T)/self.T**2

    def dQdT(self):
        """See equation 6b"""
        self.dQdT_infrared()
        self.dQdT_latentheat()
        self.dQdT_sensibleheat()
        return self.dQdT_lh + self.dQdT_ir + self.dQdT_sh
