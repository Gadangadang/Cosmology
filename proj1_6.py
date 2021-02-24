import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
from progress.bar import Bar

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


redshift = np.asarray(redshift)
lumdist = np.asarray(lumdist)*3.085678e+25#meter
lumdist_err = np.asarray(lumdist_err)*3.085678e+25#meter

w = -1
ci = c
a_0 = 1
h02 = 74
k = 1

def H2(omega_mo, omega_wo, a):
    omega_ko = 1 - omega_mo - omega_wo

    h2 = h02*omega_mo*(a_0/a)**3 + omega_ko*(a_0/a)**2 + omega_wo*(a_0/a)**(3*(1+w))
    return h2

def a(z):
    return a_0/(1+z)

def bigbang_test(z,omega_mo, omega_wo):
    H2h_02 = H2(omega_mo, omega_wo, z)
    return H2h_02 > 0

def S(k,i):
    if k == 1:
        return np.sin(i)
    elif k == 0:
        return i
    elif k == -1:
        return np.sinh(i)


def r(k,a,a_arr,H):


    i = ci*np.trapz(1/(a**2 * H ), a_arr)

    return S( k, i )

def lum_model(z,omega_mo, omega_wo):
    global k
    a_arr = np.linspace(a(z), a_0, 100)
    dl = a_0*( 1+z ) * r( k,a(z),a_arr, np.sqrt( H2(omega_mo, omega_wo, a_arr ) ) )
    return dl

def xi2(z,omega_mo, omega_wo):
    xi2 = 0
    for i in range(len(redshift)):
        #print(lumdist[i],lum_model(z,omega_mo, omega_wo),lumdist_err[i])
        xi2 += ( lumdist[i] - lum_model(z,omega_mo, omega_wo) )**2/lumdist_err[i]**2

    return xi2



z = np.linspace(np.min(redshift), np.max(redshift), 100)
omega_mo = np.linspace(0,2,100)
omega_wo = np.linspace(0,2,100)


xi2_vals = []


#Run simulation

def run():
    bar = Bar('Processing', max=len(z))
    for iz in z:

        for om_m in omega_mo:
            for om_w in omega_wo:
                xi2_vals.append( xi2(iz, om_m, om_w) )
                #if bool(bigbang_test( iz, om_m, om_w )):

        bar.next()
    bar.finish()
run()

xi2_vals = np.asarray(xi2_vals)
xi2_min = np.min(xi2_vals)

x95 = np.where(xi2_vals - xi2_min < 6.17)
print(x95)
