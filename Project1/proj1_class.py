import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
from progress.bar import Bar
import scipy.integrate as integrate
import sys


class CosmoProject1:
    def __init__(self, N=200):

        self.inter = N
        self.factor1 = 3.085678*1e25
        self.factor2 = 3.085678*1e22

        self.omega_m_min = 0
        self.omega_m_max = 1.5
        self.omega_w_min = -2
        self.omega_w_max = 3
        self.w_min = -1.5
        self.w_max = -0.5


        self.ci = c
        self.h0 = (70*1000/self.factor2)#1/s^2

    def read_file(self, fileName):
        #Read from file
        self.fileName = fileName

        file = open(self.fileName,"r")

        redshift = []
        lumdist = []
        lumdist_err = []

        for line in file:

            if len(line.split()) < 4 and len(line.split()) > 1:
                reds, lumdi, err = line.split()
                redshift.append( float( reds ) )
                lumdist.append( float( lumdi ) )
                lumdist_err.append( float( err ) )

            else:
                file.readline()



        self.redshift = np.asarray(redshift)
        self.lumdist = np.asarray(lumdist)*self.factor1  #meter
        self.lumdist_err = np.asarray(lumdist_err)*self.factor1  #meter

        file.close()


    def bigbang_test(self):
        """
        Check if H^2 < 0

        Returns
        -------
        True if H^2 > 0
        False if H^2 < 0
        """

        Hh_0 = self.h0*np.sqrt(self.current_omega_m \
                               * (1+self.redshift)**3\
                               + self.current_omega_k \
                               * (1+self.redshift)**2\
                               + self.current_omega_w \
                               * (1+self.redshift)**(3*(1+self.current_w)))**2

        if Hh_0.all() > 0:
            return True
        else:
            return False


    def S(self, k,i):
        """
        Calculate S_k(i)

        Parameters
        -----------

        k : int
          Spatial curvature parameter
        i : {array-like}, shape = [redshift]

        Returns
        -------
        sin(i), for k = 1
        i, for k = 0
        sinh(i), for k = -1
        """
        if k == 1:
            return np.sin(i)
        elif k == 0:
            return i
        elif k == -1:
            return np.sinh(i)

    def integral(self, H, z):
        """
        Calculate integral from 0 to z

        Parameters
        -----------
        H : func
          function for calculating expansion rate
        z : {array-like}, shape = [z_samples]
          Redshift values, read from file

        Returns
        -------
        res : integral values
        """
        # Do the integral for all input z
        Nz  = len(z) # Length of z
        res = np.zeros(Nz)
        res_prev = 0
        # Make sure to do the first step correctly, i.e. if z[0] is not zero
        # then integrate from 0 to z[0] first
        if(z[0] > 1e-3):
            res[0] = integrate.quad(H, 0, z[0])[0]

        for i in range(1,Nz):
            z_start  = z[i-1]
            z_end    = z[i]
            res_prev = res[i-1]
            res_next = res_prev + integrate.quad(H, z_start, z_end)[0]
            res[i]   = res_next

        return res

    def H(self, z):
        """
        Calculate 1/H as a function of z, used in r calculation

        Returns
        -------
        1/H
        """

        h = self.h0*np.sqrt(self.current_omega_m * (1+z)**3\
                            + self.current_omega_k * (1+z)**2\
                            + self.current_omega_w * (1+z)**(3 * (1+self.current_w)))
        return 1/h


    def r(self):
        """
        Calculate r integral

        Returns
        -------
        r values

        """


        i = np.sqrt(np.abs(self.current_omega_k))*self.h0

        if self.current_omega_k < 0:
            k = 1
            return self.S(k, self.integral(self.H, self.redshift)*i)
        elif self.current_omega_k > 0:
            k = -1
            return self.S(k, self.integral(self.H, self.redshift)*i)
        else:
            k = 0
            return self.integral(self.H, self.redshift)*self.h0




    def x_i2(self):
        """
        Calculate chi^2

        Returns
        -------
        xi2

        """
        #Find the numerical model for luminosity distance




        if self.current_omega_k == 0:
            lumdist_model = self.ci * ( 1+self.redshift )/(self.h0)\
                              * self.r()
        else:
            h0_sqrt_omega_k = self.h0*np.sqrt(np.abs(self.current_omega_k))

            lumdist_model = self.ci * ( 1+self.redshift )\
                              / h0_sqrt_omega_k * self.r()


        xi2 = np.sum(( self.lumdist - lumdist_model )**2/self.lumdist_err**2)
        return xi2

    def find_optimal_parameter_omega_w_omega_m_combination(self):
        """
        Finds the optimal combination of omega_w and omega_m by plotting the
        likelyhood that the combination is the true combination

        """

        self.current_w = -1
        omega_mo = np.linspace(self.omega_m_min, self.omega_m_max, self.inter)
        omega_wo = np.linspace(self.omega_w_min, self.omega_w_max, self.inter)

        X,Y = np.meshgrid(omega_mo, omega_wo)
        self.xi2val = np.zeros( ( self.inter, self.inter ) )


        #Run simulation
        #Remove hashtag in the run function as well as in
        #the package section to have progress bar in terminal

        #bar = Bar('Processing', max=self.inter)

        for iom,om_m in enumerate(omega_mo):
            self.current_omega_m = om_m

            for iow,om_w in enumerate(omega_wo):

                self.current_omega_w = om_w
                self.current_omega_k = 1 - self.current_omega_m\
                                         - self.current_omega_w

                """
                print("o_m = ", self.current_omega_m)
                print("o_w = ", self.current_omega_w)
                print("o_k = ", self.current_omega_k)"""


                if self.bigbang_test() is True:
                    val = self.x_i2()
                    self.xi2val[iow, iom] = val
                else:
                    self.xi2val[iow, iom] = np.nan

            #bar.next()
        #bar.finish()

        self.xi2val = self.xi2val - np.nanmin(self.xi2val)

        cmap = plt.cm.coolwarm
        cmap.set_bad('black',1.)
        plt.contourf(X,Y,np.log10(self.xi2val),cmap=cmap)

        plt.colorbar()
        plt.xlabel(r"$\Omega_{m 0}$")
        plt.ylabel(r"$\Omega_{w 0}$")
        plt.savefig("images/totxi2_6.jpeg")
        plt.show()



        # 95% range
        index = np.where(self.xi2val - np.nanmin(self.xi2val) < 6.17,\
                         self.xi2val, np.nan)

        cmap = plt.cm.coolwarm
        cmap.set_bad('black',1.)
        plt.contourf(X,Y,np.log10(index),cmap=cmap)
        plt.tight_layout()
        plt.colorbar()

        plt.xlabel(r"$\Omega_{m 0}$")
        plt.ylabel(r"$\Omega_{w 0}$")

        plt.savefig("images/95%xi2_6.jpeg")
        plt.show()

    def find_optimal_parameter_w_omega_w_combination(self):
        """
        Finds the optimal combination of omega_w and omega_m by plotting the
        likelyhood that the combination is the true combination

        """

        self.current_omega_m = 0.3
        w = np.linspace(self.w_min, self.w_max, self.inter)
        omega_wo = np.linspace(self.omega_w_min, self.omega_w_max, self.inter)

        X,Y = np.meshgrid(omega_wo, w)
        self.xi2val = np.zeros( ( self.inter, self.inter ) )


        #Run simulation
        #Remove hashtag in the run function as well as in
        #the package section to have progress bar in terminal

        #bar = Bar('Processing', max=self.inter)

        for iw,w_val in enumerate(w):
            self.current_w = w_val

            for iow,om_w in enumerate(omega_wo):

                self.current_omega_w = om_w
                self.current_omega_k = 1 - self.current_omega_m\
                                         - self.current_omega_w

                """
                print("o_m = ", self.current_omega_m)
                print("o_w = ", self.current_omega_w)
                print("o_k = ", self.current_omega_k)"""


                if self.bigbang_test() is True:
                    val = self.x_i2()
                    self.xi2val[iw, iow] = val
                else:
                    self.xi2val[iw, iow] = np.nan



        self.xi2val = self.xi2val - np.nanmin(self.xi2val)

        cmap = plt.cm.coolwarm
        cmap.set_bad('black',1.)
        plt.contourf(X,Y,np.log10(self.xi2val),cmap=cmap)

        plt.colorbar()
        plt.xlabel(r"$\Omega_{m 0}$")
        plt.ylabel(r"$\Omega_{w 0}$")
        plt.savefig("images/totxi2_7.jpeg")
        plt.show()



        # 95% range
        index = np.where(self.xi2val - np.nanmin(self.xi2val) < 6.17,\
                         self.xi2val, np.nan)

        cmap = plt.cm.coolwarm
        cmap.set_bad('black',1.)
        plt.contourf(X,Y,np.log10(index),cmap=cmap)
        plt.tight_layout()
        plt.colorbar()

        plt.xlabel(r"$\Omega_{m 0}$")
        plt.ylabel(r"$\Omega_{w 0}$")

        plt.savefig("images/95%xi2_7.jpeg")
        plt.show()




if "__main__" == __name__:
    Omega = CosmoProject1(N=40)
    Omega.read_file("sndata.txt")
    Omega.find_optimal_parameter_omega_w_omega_m_combination()
    Omega.find_optimal_parameter_w_omega_w_combination()
