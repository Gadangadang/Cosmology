import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as intp
from scipy import constants
import scipy.misc



class TemperatureTable:

    def __init__(self, T):
        """
        Initialize class variables and
        arrays used in calculation
        """

        self.Temp = T
        self.m = constants.electron_mass
        self.hbar = constants.hbar
        self.c = constants.c
        self.kb = constants.k
        self.G = constants.G

        self.T_neutrino = np.zeros_like(self.Temp)
        self.t = np.zeros_like(self.Temp)
        self.Splot = np.zeros_like(self.Temp)

    def inside_integral_S(self, y, x):
        """
        Algebraic expression inside integral in S(x)

        Parameters
        -----------

        y : float
          Integral parameter
        x : float
          constant given as mc^2/kT

        Returns
        -------
        Expression
        """
        return y**2 * (np.sqrt( x**2 + y**2 ) + y**2 \
               / ( 3 * np.sqrt( x**2 + y**2 ))) \
               / ( np.exp( np.sqrt( x**2 + y**2 )) + 1 )

    def inside_integral_E(self, y, x):
        """
        Algebraic expression inside integral in E(x)

        Parameters
        -----------

        y : float
          Integral parameter
        x : float
          constant given as mc^2/kT

        Returns
        -------
        Expression
        """
        return y**2 * (np.sqrt(x**2 + y**2)) \
               / (np.exp( np.sqrt(x**2 + y**2) ) + 1)


    def find_epsilon(self, x):
        """
        Calculate E(x)

        Parameters
        -----------

        x : float
          constant given as mc^2/kT

        Returns
        -------
        Solved integral
        """
        epint = 30/np.pi**4 * float(intp.quad(self.inside_integral_E,
                                              0,
                                              np.inf,
                                              args=(x,))[0])

        epsi = 1 + 21/8 * (4/11)**(4/3) * self.find_S_of_x(x)**(4/3) + epint
        return epsi

    def inside_integral_t(self, x):
        """
        Algebraic expression inside integral in t(s)

        Parameters
        -----------

        x : float
          constant given as mc^2/kT

        Returns
        -------
        Expression
        """
        return (3 - x*scipy.misc.derivative(self.find_S_of_x, x)\
                /self.find_S_of_x(x))*self.find_epsilon(x)**(-1/2)*x

    def find_S_of_x(self, x):
        """
        Calculate S(x)

        Parameters
        -----------

        x : float
          constant given as mc^2/kT

        Returns
        -------
        Solved S(x) value
        """
        inter = float(intp.quad(self.inside_integral_S,
                                0,
                                np.inf,
                                args=(x,))[0])

        return 1 + inter*45/(2*np.pi**4)

    def find_x_of_T(self, T):
        """
        Calculate x(T)

        Parameters
        -----------

        T : float

        Returns
        -------
        x(T)
        """
        return self.m*self.c**2/(self.kb*T)

    def find_t_of_s(self, T):
        """
        Calculate t(s)

        Parameters
        -----------

        T : float
          Temperature [K]

        Returns
        -------
        Array of time values
        """
        x_vals = self.find_x_of_T(T)
        x0 = x_vals[0]
        const = np.sqrt((15 * self.hbar**3)\
                         / (24 * np.pi**3 * self.G \
                         * self.m**4 * self.c**3) )

        xlength = len(x_vals)
        results = np.zeros(xlength)

        for i,xi in enumerate(x_vals):
            results[i] = intp.quad(self.inside_integral_t, x0, xi)[0]

        return const*results

    def find_neutrino_temperature(self, T):
        """
        Calculate T_neutrino

        Parameters
        -----------

        T : float
          Temperature T [K]

        Returns
        -------
        T_neutrino value
        """
        T_neu_val = (4/11)**(1/3)*T\
                    *self.find_S_of_x(self.find_x_of_T(T))**(1/3)
        return T_neu_val

    def calc_and_print(self):
        """
        Calculate T_neutrino, t(s) and print table
        """

        self.t = self.find_t_of_s(self.Temp)

        print("                                     ")
        print("     T       T_n/T       t      ")
        print("-------------------------------------")

        for i in range(len(self.Temp-1)):
            self.T_neutrino[i] = self.find_neutrino_temperature(self.Temp[i])
            self.Splot[i] = self.find_S_of_x(self.find_x_of_T(self.Temp[i]))

            print("{:.2e}  |  {:.5f}  |  {:.4e} "\
                  .format(self.Temp[i],
                          self.T_neutrino[i]/self.Temp[i],
                          self.t[i]))

    def plots(self):
        """
        Plot S(x) and T_neutrino/T (T)
        """
        plt.scatter(np.log10(self.find_x_of_T(self.Temp)),
                    self.Splot,
                    label="S(x)")
        plt.legend()
        plt.xlabel(r"$log_{10}(X=m_ec^2/k_BT)$")
        plt.ylabel(r"$S(X=m_ec^2/k_BT)$")
        plt.savefig("Images/s_x.jpeg")
        plt.show()

        plt.plot(self.Temp,
                 self.T_neutrino/self.Temp,
                 label=r"$T_{\nu}/T$")
        plt.legend()
        plt.xlabel("T[K]")
        plt.ylabel(r"$T_{\nu}/T$")
        plt.savefig("Images/T_T.jpeg")
        plt.show()



if __name__ == "__main__":

    #Test temperatures
    T = np.array([1e11, 6e10, 2e10,
                  1e10, 6e9, 3e9,
                  2e9, 1e9, 3e8,
                  1e8, 1e7, 1e6])

    calcis = TemperatureTable(T)
    anal_x0 = 1 + 45/(2*np.pi**4)*4/3*7*np.pi**4/120

    #Check model for x = 0 and x >> 1
    print("For x = 0, model gives {:.3f} and analytical gives {:.3f}"\
          .format(calcis.find_S_of_x(0), anal_x0))
    print("For x >> 1, model gives {:.3f} and analytical gives {:.3f}"\
          .format(calcis.find_S_of_x(50), 1.))

    #Find values for T_neutrino and t(s)
    calcis.calc_and_print()
    calcis.plots()
