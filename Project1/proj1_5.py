import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
import scipy.integrate as integrate
import sys

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
    omega_ko = 1 - omega_m - omega_w

    if omega_ko == 0:
        dl = c*( 1+z )/(h0) * r( z, omega_ko, omega_m, omega_w)
    else:
        dl = c*( 1+z )/(h0*np.sqrt(np.abs(omega_ko))) * r( z, omega_ko, omega_m, omega_w)

    return dl

factor2 = 3.085678*1e22
h0 = (70)#1/s^2
z = np.linspace(0, 2, 1001)
c = c



#Example 1 Einstein - de Sitter model
def Eds_ana(z):
    return 2*c/h0 * (1 + z - np.sqrt(1 + z))


w = 0
omega_m = 1
omega_w = 0

numdl_eds = lum_model(z, omega_m, omega_w)

#Example 2 Milne model
def Milne_ana(z):
    return c*z/(2*h0) *(2 + z)


w = -1/3
omega_m0 = 0
omega_w0 = 0



numdl_milne = lum_model(z, omega_m0, omega_w0)

#Plot solutions

plt.plot(z, Eds_ana(z)/1e6, label="Analytical E-dS")
plt.plot(z, numdl_eds/1e6,"--",label="Numerical E-dS")
plt.plot(z, Milne_ana(z)/1e6, label="Analytical Milne")
plt.plot(z, numdl_milne/1e6, "--", label="Numerical Milne")

plt.ylabel(r"$d_L(z)$ [MPc]")
plt.xlabel("z redshift")
plt.legend()
plt.savefig("images/5_testmodel.jpeg")
plt.show()

#Relative error as a function of z

relerr_EdS = (Eds_ana(z) - numdl_eds)/Eds_ana(z)
relerr_Milne = (Milne_ana(z) - numdl_milne)/Milne_ana(z)

plt.plot(z, relerr_EdS, label="Relative Error E - dS")
plt.plot(z, relerr_Milne, label="Relative Error Milne")

plt.ylabel("Relative error")
plt.xlabel("z redshift")
plt.legend()
plt.savefig("images/5_relerror.jpeg")
plt.show()
