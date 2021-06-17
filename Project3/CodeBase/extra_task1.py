import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

c = constants.c
G = constants.G

solar_mass = 5.0289921396853E-31

mass_a = 170 / solar_mass
mass_b = 170 / solar_mass

N = 2/5 #EM correction
eta = (mass_a*mass_b)/(mass_a+mass_b)**2

R_distance_to_us = 440*3.0856776e+22 #Mpc to meter
R_s = 2*G*(mass_a + mass_b)/c**2
init_r_separation = 1.7*R_s




def r(t):
    return (((init_r_separation)**4 - (N*eta*c*t*R_s**3))**(1/4))

def omega(t):
    return (c/np.sqrt(2))*np.sqrt( R_s/r(t)**3 )


def h(t):
    return ((1/2)*(1/R_distance_to_us)*mass_a*mass_b*G**2/(2*c**4*r(t)))\
            *np.cos(2*omega(t)*t)


T = 0.27
N_w = 10000
dt = T/N_w
strain_values = np.zeros(N_w)
time = np.linspace(0.0,0.27,N_w)
i = 0

while time[i] < T:
    strain_values[i] = h(time[i])
    i = i + 1



strain_values = strain_values/np.max(strain_values) #Normalize strain values

plt.plot(time, strain_values)
plt.ylabel("h(t)")
plt.xlabel("t [s]")
plt.title(r"$\frac{h(t)}{h_{max}}$ strain GW170104 observation")
plt.savefig("figure6_remake.jpeg")
plt.show()
