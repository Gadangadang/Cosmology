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
        self.G = constants.G

        self.T_neutrino = np.zeros_like(self.Temp)
        self.t = np.zeros_like(self.Temp)
        self.Splot = np.zeros_like(self.Temp)

    def inteS_x(self, y, x):
        return y**2 * (np.sqrt( x**2 + y**2 ) + y**2 \
               / ( 3 * np.sqrt( x**2 + y**2 ))) \
               / ( np.exp( np.sqrt( x**2 + y**2 )) + 1 )

    def inteE_x(self, y, x):
        return y**2 * (np.sqrt(x**2 + y**2)) \
               / (np.exp( np.sqrt(x**2 + y**2) ) + 1)

    def S_der_x_int(self, y, x):
        return -(y**2 * x * ((3 * (x**2 + y**2)**(3/2) + y**2 \
               * np.sqrt(x**2 + y**2) - 3 * x**2 - 2 * y**2) \
               * np.exp(np.sqrt(x**2 + y**2)) - 3 * x**2 - 2 * y**2)) \
               / (3 * (x**2+y**2)**(3/2) \
               * (np.exp(np.sqrt(x**2 + y**2)) + 1)**2)

    def S_der(self, x):
        inter = float(intp.quad(self.S_der_x_int, 0, np.inf, args=(x))[0])

        return inter*45/(2*np.pi**4)

    def epsilon(self, x):
        epint = 30/np.pi**4 * float(intp.quad(self.inteE_x,
                                              0,
                                              np.inf,
                                              args=(x))[0])

        epsi = 1 + 21/8 * (4/11)**(4/3) * self.S(x)**(4/3) + epint
        return epsi

    def t_int_x(self, x):
        return (3 - x*self.S_der(x)/self.S(x))*self.epsilon(x)**(-1/2)*x

    def S(self, x):
        inter = float(intp.quad(self.inteS_x, 0, np.inf, args=(x))[0])
        return 1 + inter*45/(2*np.pi**4)

    def x(self, T):
        return self.m*self.c**2/(self.kb*T)

    def t_vals(self, T):
        x_vals = self.x(T)


        const = np.sqrt((15 * self.hbar**3)\
                        / (24 * np.pi**3 * self.G \
                        * self.m**4 * self.c**3) )



        xlength = len(x_vals)
        results = np.zeros(xlength)
        res_prev = 0


        for i in range(xlength):
            x_start  = x_vals[i-1]
            x_end    = x_vals[i]
            res_prev = results[i-1]
            res_next = res_prev + intp.quad(self.t_int_x, x_start, x_end)[0]
            results[i] = res_next


        return const*results

    def T_neut(self, T):
        T_neu_val = (4/11)**(1/3)*T*self.S(self.x(T))**(1/3)
        return T_neu_val

    def calc_and_print(self):

        self.t = self.t_vals(self.Temp)

        print("     T       T_neutrino       t      ")
        print("-------------------------------------")

        for i in range(len(self.Temp-1)):
            self.T_neutrino[i] = self.T_neut(self.Temp[i])
            self.Splot[i] = self.S(self.x(self.Temp[i]))

            print("{:.3e}  |  {:.3e}  |  {:.3e} "\
                  .format(self.Temp[i],
                          self.T_neutrino[i]/self.Temp[i],
                          self.t[i]))

    def plots(self):

        plt.scatter(self.x(self.Temp), self.Splot, label="S(x)")
        plt.legend()
        plt.savefig("Images/s_x.jpeg")
        plt.show()



if __name__ == "__main__":
    T = np.array([1e11, 6e10, 2e10, 1e10, 6e9, 3e9, 2e9, 1e9, 3e8, 1e8, 1e7, 1e6])
    calcis = TemperatureTable(T)

    print(calcis.S(0), 1 + 45/(2*np.pi**4)*7*np.pi**4/120)
    calcis.calc_and_print()
    calcis.plots()
