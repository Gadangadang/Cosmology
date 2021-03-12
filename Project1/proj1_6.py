import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
from progress.bar import Bar
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
a_0 = 1
h0 = (70*1000/factor2)#1/s^2





#Functions for analysis



def a(z):
    return a_0/(1+z)

def bigbang_test(z,omega_m, omega_w):
    omega_ko = 1 - omega_m - omega_w
    Hh_0 = h0*np.sqrt(omega_m*(1+z)**3 + omega_ko*(1+z)**2 + omega_w*(1+z)**(3*(1+w)))**2


    if Hh_0.all() > 0:
        return True
    else:
        return False

def S(k,i):
    if k == 1:
        return np.sin(i)
    elif k == 0:
        return i
    elif k == -1:
        return np.sinh(i)





def r(redshift, omega_ko, omega_m, omega_w):

    def integral(H, z):


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
    global redshift
    global size
    omega_ko = 1 - omega_m - omega_w

    if omega_ko == 0:
        dl = c*( 1+z )/(h0) * r( redshift, omega_ko, omega_m, omega_w  )
    else:
        dl = c*( 1+z )/(h0*np.sqrt(np.abs(omega_ko))) * r( redshift, omega_ko, omega_m, omega_w  )

    return dl

def xi2(z,omega_m, omega_w):

    xi2 = np.sum( ( lumdist - lum_model(z,omega_m, omega_w) )**2/lumdist_err**2 )
    return xi2

size = len(redshift)

inter = 300

omega_mo = np.linspace(0,1.5,inter)
omega_wo = np.linspace(0,1.5,inter)

X,Y = np.meshgrid(omega_mo, omega_wo)
xi2val = np.zeros( ( len(omega_mo),len(omega_wo) ) )



#Run simulation

def run():
    bar = Bar('Processing', max=inter)

    for iom,om_m in enumerate(omega_mo):
        for iow,om_w in enumerate(omega_wo):
            if bigbang_test(redshift, om_m, om_w) is True:
                val = xi2(redshift, om_m, om_w)
                xi2val[iow, iom] = val
            else:
                xi2val[iow, iom] = np.nan

        bar.next()
    bar.finish()

run()

xi2val = xi2val - np.nanmin(xi2val)
#rint(xi2val)

#print(np.min(xi2val))

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
