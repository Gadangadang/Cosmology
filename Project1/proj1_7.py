import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
#from progress.bar import Bar
import scipy.integrate as integrate

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

def bigbang_test(z,omega_w, w):
    """
    Check if H^2 < 0

    Parameters
    -----------

    z : {array-like}, shape = [z_samples]
      Redshift values, read from file
    w : float
      w value
    omega_w : float
      Omega_w0 value

    Returns
    -------
    True if H^2 > 0
    False if H^2 < 0
    """
    global omega_m
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





def r(redshift, omega_ko, w, omega_w):
    """
    Calculate r integral

    Parameters
    -----------
    z : {array-like}, shape = [z_samples]
      Redshift values, read from file
    omega_k0 : float
      omega_k0 value to calculate r integral
    w : float
      w value to calculate r integral
    omega_w : float
      omega_w value to calculate integral

    Returns
    -------
    r values

    """
    global omega_m

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





def lum_model(z,w, omega_w):
    """
    Calculate luminosity distance

    Parameters
    -----------
    z : {array-like}, shape = [z_samples]
      Redshift values, read from file
    w : float
      w value to calculate luminosity distance
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
        dl = c*( 1+z )/(h0*np.sqrt(np.abs(omega_ko))) * r( redshift, omega_ko, w, omega_w  )

    return dl

def xi2(z,w, omega_w):
    """
    Calculate chi^2

    Parameters
    -----------
    z : {array-like}, shape = [z_samples]
      Redshift values, read from file
    w : float
      w value to calculate xi2
    omega_w : float
      omega_w value to calculate xi2

    Returns
    -------
    xi2

    """

    xi2 = np.sum( ( lumdist - lum_model(z,w, omega_w) )**2/lumdist_err**2 )
    return xi2

# Set grid for computing a contourf plot
size = len(redshift)
inter = 300
omega_m = 0.3
omega_wo = np.linspace(-2,3,inter)
w = np.linspace(-1.5,-0.5,inter)

X,Y = np.meshgrid(omega_wo, w)
xi2val = np.zeros( ( len(omega_wo),len(w) ) )

#Run simulation
#Remove hashtag in the run function as well as in
#the package section to have progress bar in terminal

def run():
    #bar = Bar('Processing', max=inter)

    for iw,w_val in enumerate(w):
        for iow,om_w in enumerate(omega_wo):
            if bigbang_test(redshift, om_w, w_val) is True:
                val = xi2(redshift, w_val, om_w)
                xi2val[iw, iow] = val


            else:
                print("p")
                xi2val[iow, iw] = np.nan


        #bar.next()
    #bar.finish()

run()


#Plot

xi2val = xi2val - np.nanmin(xi2val)


cmap = plt.cm.coolwarm
cmap.set_bad('black',1.)
plt.contourf(X,Y,np.log10(xi2val),cmap=cmap)

plt.colorbar()
plt.xlabel(r"$\Omega_{w 0}$")
plt.ylabel("w")
plt.savefig("images/xi2_7.jpeg")
plt.show()

file.close()

# 95%

index = np.where( xi2val - np.nanmin(xi2val) < 6.17, xi2val, np.nan )
#rint(xi2val)

#print(np.min(xi2val))



cmap = plt.cm.coolwarm
cmap.set_bad('black',1.)
plt.contourf(X,Y,np.log10(index),cmap=cmap)

plt.colorbar()
plt.xlabel(r"$\Omega_{w 0}$")
plt.ylabel("w")

plt.savefig("images/95%xi2_7.jpeg")
plt.show()

file.close()
