import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

G = 4*constants.pi**2 #solar masses
c = constants.c*0.00021 #AY/yr

def schwartzchild_radius(M):
    return 2*G*M/c**2

M = 70 #solar masses
r_i = 4.7*schwartzchild_radius(M) # AU and solar mass units
t_i = 0
eta = 1/4 # two equal solar masses
N = 32/5 # GR model

yr = 3.16e-8 #sec to yr

t = np.linspace(t_i, 0.2, 1001)

def radius_separation_over_rs(t):
    term = (( (r_i)/(schwartzchild_radius(M)))**(4) \
              - (N*eta*c*(t-t_i))/(schwartzchild_radius(M)) )**(1/4)
    return term

def omega(t):
    term = (c/(np.sqrt(2)))*np.sqrt( schwartzchild_radius(M) \
                                   /(radius_separation_over_rs(t)\
                                   *schwartzchild_radius(M))**3 )
    return term

def v(t):
    return radius_separation_over_rs(t)*omega(t)*schwartzchild_radius(M)/2


fig,ax = plt.subplots()

ax.plot(t, radius_separation_over_rs(t*yr), color="Black")
ax.set_ylabel("Separation", color="black")
ax.set_xlabel("t[s]")
ax2=ax.twinx()
ax2.plot(t, v(t*yr)/c, color="Blue")
ax2.set_ylabel("Speed/c", color="Blue")

fig.savefig("separation_and_speed_for_35Msol_black_holes.jpg",
            format='jpeg',
            dpi=100,
            bbox_inches='tight')

plt.show()

ma = 20
mb = 35
M = ma + mb
eta = (ma*mb)/(ma+mb)**2
print(eta)


fig,ax = plt.subplots()

ax.plot(t, radius_separation_over_rs(t*yr), color="Black")
ax.set_ylabel("Separation", color="black")
ax.set_xlabel("t[s]")
ax2=ax.twinx()
ax2.plot(t, v(t*yr)/c, color="Blue")
ax2.set_ylabel("Speed/c", color="Blue")

fig.savefig("separation_and_speed_for_20_35_Msol_black_holes.jpg",
            format='jpeg',
            dpi=100,
            bbox_inches='tight')

plt.show()

ma = 1
mb = 30
M = ma + mb
eta = (ma*mb)/(ma+mb)**2
print(eta)

fig,ax = plt.subplots()

ax.plot(t, radius_separation_over_rs(t*yr), color="Black")
ax.set_ylabel("Separation", color="black")
ax.set_xlabel("t[s]")
ax2=ax.twinx()
ax2.plot(t, v(t*yr)/c, color="Blue")
ax2.set_ylabel("Speed/c", color="Blue")

fig.savefig("separation_and_speed_for_30_1_Msol_black_holes.jpg",
            format='jpeg',
            dpi=100,
            bbox_inches='tight')

plt.show()
