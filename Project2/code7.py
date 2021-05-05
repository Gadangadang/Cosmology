import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as intp
from scipy import constants



class TemperatureTable:
    def __init__(self, T):
        self.Temp = T
        self.m = constants.electron_mass
        self.hbar = constants.hbar
        self.c = constants.c
        self.kb = constants.k

    def inteS_x(self, y, x):
        return y**2 * (np.sqrt( x**2 + y**2 )\
               + y**2 / ( 3 * np.sqrt( x**2 + y**2 )))\
               / ( np.exp( np.sqrt( x**2 + y**2 )) + 1 )

    def S(self,x):
        inter = 45/(2*np.pi**4)
        inter *= float(intp.quad(self.inteS_x, 0, np.inf, args=(x))[0])
        return 1 + inter

    def x(self,T):
        return self.m*self.hbar/(self.kb*T)

    def t_vals(self, T):
        return 0

    def T_neut(self, T):
        T_neu_val = (4/11)**(1/3)*T*self.S(self.x(T))**(1/3)
        return T_neu_val

    def calc_and_print(self):
        self.T_neutrino = np.zeros_like(self.Temp)
        self.t = np.zeros_like(self.Temp)

        print("     T       T_neutrino       t      ")
        print("-------------------------------------")

        for i in range(len(self.Temp-1)):
            self.T_neutrino[i] = self.T_neut(self.Temp[i])
            self.t[i] = self.t_vals(self.Temp[i])

            print("{:.3e}  |  {:.3e}  |  {:.3e} ".format(self.Temp[i],
                                                         self.T_neutrino[i],
                                                         self.t[i]))

if __name__ == "__main__":
    T = np.array([1e11, 6e10, 2e10, 1e10, 6e9, 3e9, 2e9, 1e9, 3e8, 1e8, 1e7, 1e6])
    calcis = TemperatureTable(T)
    print(calcis.S(0), 1 + 45/(2*np.pi**4)*7*np.pi**4/120)
    calcis.calc_and_print()
