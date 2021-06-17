import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

c = constants.c
G = constants.G

solar_mass = 5.0289921396853E-31

mass_a = 1.4 / solar_mass
mass_b = 1.9 / solar_mass

N = 2/5 #EM correction
eta = (mass_a*mass_b)/(mass_a+mass_b)**2

R_distance_to_us = 200/1.0570008340246E-16 #200 Ly to meters
R_s = 2*G*(mass_a + mass_b)/c**2
init_r_separation = 9*R_s




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

print("Strain strength for binary neutron system is {:.3e}".format(np.max(strain_values)))
print("-------------")
strain_values = strain_values/np.max(strain_values) #Normalize strain values



plt.plot(time, strain_values)
plt.ylabel("h(t)")
plt.xlabel("t [s]")
plt.title(r"$\frac{h(t)}{h_{max}}$ strain GW170104 observation")
plt.savefig("h_strain_of_neutron star.jpeg")
plt.show()

plt.subplot(211)
plt.title(r"$\omega(t)$ and $\dot{\omega}(t)$ for binary neutron stars")
plt.plot(time,omega(time),label=r"$\omega(t)$")
plt.xlabel("t [s]")
plt.ylabel(r"$\omega(t)$ [Hz]")

plt.legend()
plt.subplot(212)

plt.plot(time, 1/2*omega(time)*r(time)/c,label=r"$ v(t)/c $")
plt.xlabel("t [s]")
plt.ylabel(r"$v(t)/c$ [m/s]")


plt.legend()
plt.savefig("omega_and_v.jpeg")
plt.show()



while np.max(strain_values) > 1e-22:
    R_distance_to_us = R_distance_to_us * 2
    strain_values = np.zeros(N_w)

    i = 0
    time = np.linspace(0,T,N_w)
    while time[i] < T:


        strain_values[i] = h(time[i])
        i = i + 1

    if np.max(strain_values) > 1e-22:
        if R_distance_to_us*1.0570008340246E-16*3.0660139380653215e-7 > 0.05:
            print("Distance is {:.3e} Light years, or {:.3f} Mpc".format(R_distance_to_us*1.0570008340246E-16, R_distance_to_us*1.0570008340246E-16*3.0660139380653215e-7))
            print("Strain strength is {:.3e}".format(np.max(strain_values)))
        else:
            print("Distance is {:.3e} Light years, or {:.3f} ly".format(R_distance_to_us*1.0570008340246E-16, R_distance_to_us*1.0570008340246E-16))
            print("Strain strength is {:.3e}".format(np.max(strain_values)))
    else:
        pass
print("-------------")
