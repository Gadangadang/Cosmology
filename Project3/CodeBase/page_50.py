import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from scipy.misc import derivative


c = constants.c
G = constants.G
#GW170104
solar_mass = 5.0289921396853E-31
mass_a = 30.8 / solar_mass #kg
mass_b = 20.0 / solar_mass #kg
N = 2/5 #EM correction
eta = (mass_a*mass_b)/(mass_a+mass_b)**2

R_distance_to_us = 3.054820824e+25 #990Mpc
R_s = 2*G*(mass_a + mass_b)/c**2
init_r_separation = 3.5*R_s




def r(t):
    return (((init_r_separation)**4 - (N*eta*c*t*R_s**3))**(1/4))

def omega(t):
    return (c/np.sqrt(2))*np.sqrt( R_s/r(t)**3 )


def h(t):
    return ((1/2)*(1/R_distance_to_us)*mass_a*mass_b*G**2/(2*c**4*r(t)))\
            *np.cos(2*omega(t)*t)


T = 0.5
N_w = 10000
dt = T/N_w
strain_values = np.zeros(N_w)
time = np.linspace(0,T,N_w)
i = 0

while time[i] < T:
    strain_values[i] = h(time[i])
    i = i + 1

omeg = omega(time)
omega_dot = np.diff(omeg) #Omega dot values

strain_values = strain_values/np.max(strain_values) #Normalize strain values

plt.plot(time, strain_values)
plt.ylabel("h(t)")
plt.xlabel("t [s]")
plt.title(r"$\frac{h(t)}{h_{max}}$ strain GW170104 observation")
plt.savefig("h_strain_of_GW170104.jpeg")
plt.show()

plt.subplot(211)
plt.title(r"$\omega(t)$ and $\dot{\omega}(t)$ for GW170104 observation")
plt.plot(time,omega(time),label=r"$\omega(t)$")
plt.xlabel("t [s]")
plt.ylabel(r"$\omega(t)$ [Hz]")

plt.legend()
plt.subplot(212)

plt.plot(time[:-1], omega_dot,label=r"$ \dot{\omega}(t) $")
plt.xlabel("t [s]")
plt.ylabel(r"$\dot{\omega}(t)$ [Hz]")


plt.legend()
plt.savefig("omega_and_omegadot.jpeg")
plt.show()
