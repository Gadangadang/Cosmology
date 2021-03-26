import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
#from progress.bar import Bar
import scipy.integrate as integrate
import sys


#Read from file

file = open("sndata.txt","r")

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

factor1 = 3.085678*1e25
factor2 = 3.085678*1e22

redshift = np.asarray(redshift)
lumdist = np.asarray(lumdist)*factor1  #meter
lumdist_err = np.asarray(lumdist_err)*factor1  #meter



w = -1
ci = c
h0 = (70*1000/factor2)#1/s^2





#Functions for analysis


def bigbang_test(z,omega_m, omega_w):
    """
    Check if H^2 < 0

    Parameters
    -----------

    z : {array-like}, shape = [z_samples]
      Redshift values, read from file
    omega_m : float
      Omega_m0 value
    omega_w : float
      Omega_w0 value

    Returns
    -------
    True if H^2 > 0
    False if H^2 < 0
    """
    omega_ko = 1 - omega_m - omega_w
    Hh_0 = h0*np.sqrt(omega_m*(1+z)**3 + omega_ko*(1+z)**2 + omega_w*(1+z)**(3*(1+w)))**2


    if Hh_0.all() > 0:
        return True
    else:
        return False

def S(k,i):
    """
    Calculate S_k(i)

    Parameters
    -----------

    k : int
      Spatial curvature parameter
    i : {array-like}, shape = [z_samples]

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





def r(redshift, omega_ko, omega_m, omega_w):
    """
    Calculate r integral

    Parameters
    -----------
    z : {array-like}, shape = [z_samples]
      Redshift values, read from file
    omega_k0 : float
      omega_k0 value to calculate r integral
    omega_m : float
      omega_m value to calculate r integral
    omega_w : float
      omega_w value to calculate integral

    Returns
    -------
    r values

    """

    def integral(H, z):
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
        # Make sure to do the first step correctly, i.e. if z[0] is not zerom then integrate from 0 to z[0] first
        if(z[0] > 1e-3):
            res[0] = integrate.quad(H, 0, z[0])[0]

        for i in range(1,Nz):
            z_start  = z[i-1]
            z_end    = z[i]
            res_prev = res[i-1]
            res_next = res_prev + integrate.quad(H, z_start, z_end)[0]
            res[i]   = res_next

        return res

    def H(z):
        """
        Calculate 1/H as a function of z, used in r calculation

        Parameters
        -----------
        z : {array-like}, shape = [z_samples]
          Redshift values, read from file

        Returns
        -------
        1/H
        """

        h = h0*np.sqrt(omega_m*(1+z)**3 + omega_ko*(1+z)**2 + omega_w*(1+z)**(3*(1+w)))
        return 1/h

    i = np.sqrt(np.abs(omega_ko))*h0

    if omega_ko < 0:
        k = 1
        return S(k, integral(H, redshift)*i)
    elif omega_ko > 0:
        k = -1
        return S(k, integral(H, redshift)*i)
    else:
        k = 0
        return integral(H, redshift)*h0





def lum_model(z,omega_m, omega_w):
    """
    Calculate luminosity distance

    Parameters
    -----------
    z : {array-like}, shape = [z_samples]
      Redshift values, read from file
    omega_m : float
      omega_m value to calculate luminosity distance
    omega_w : float
      omega_w value to calculate luminosity distance

    Returns
    -------
    dl : luminosity distance

    """
    omega_ko = 1 - omega_m - omega_w

    if omega_ko == 0:
        dl = c*( 1+z )/(h0) * r( redshift, omega_ko, omega_m, omega_w  )
    else:
        dl = c*( 1+z )/(h0*np.sqrt(np.abs(omega_ko))) * r( redshift, omega_ko, omega_m, omega_w  )

    return dl

def xi2(z,omega_m, omega_w):
    """
    Calculate chi^2

    Parameters
    -----------
    z : {array-like}, shape = [z_samples]
      Redshift values, read from file
    omega_m : float
      omega_m value to calculate xi2
    omega_w : float
      omega_w value to calculate xi2

    Returns
    -------
    xi2

    """

    xi2 = np.sum( ( lumdist - lum_model(z,omega_m, omega_w) )**2/lumdist_err**2 )
    return xi2


# Set grid for computing a contourf plot
size = len(redshift)
inter = 200

omega_mo = np.linspace(0,1.5,inter)
omega_wo = np.linspace(-2,3,inter)

X,Y = np.meshgrid(omega_mo, omega_wo)
xi2val = np.zeros( ( len(omega_mo),len(omega_wo) ) )



#Run simulation
#Remove hashtag in the run function as well as in
#the package section to have progress bar in terminal


def run():
    #bar = Bar('Processing', max=inter)

    for iom,om_m in enumerate(omega_mo):
        for iow,om_w in enumerate(omega_wo):
            if bigbang_test(redshift, om_m, om_w) is True:
                val = xi2(redshift, om_m, om_w)
                xi2val[iow, iom] = val
            else:
                xi2val[iow, iom] = np.nan

        #bar.next()
    #bar.finish()

run()

xi2val = xi2val - np.nanmin(xi2val)


#Plot code

cmap = plt.cm.coolwarm
cmap.set_bad('black',1.)
plt.contourf(X,Y,np.log10(xi2val),cmap=cmap)

plt.colorbar()
plt.xlabel(r"$\Omega_{m 0}$")
plt.ylabel(r"$\Omega_{w 0}$")
plt.savefig("images/totxi2_6.jpeg")
plt.show()

file.close()

# 95% range
index = np.where(xi2val - np.nanmin(xi2val) < 6.17, xi2val, np.nan)


#sys.exit()

cmap = plt.cm.coolwarm
cmap.set_bad('black',1.)
plt.contourf(X,Y,np.log10(index),cmap=cmap)
plt.tight_layout()
plt.colorbar()

plt.xlabel(r"$\Omega_{m 0}$")
plt.ylabel(r"$\Omega_{w 0}$")

plt.savefig("images/95%xi2_6.jpeg")
plt.show()

file.close()
