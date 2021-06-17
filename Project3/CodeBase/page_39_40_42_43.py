from scipy import constants
import numpy as np


G = constants.G
c = constants.c
omega1 = 2*np.pi/(0.275-0.25)
omega2 = 2*np.pi/(0.298-0.275)
omega = (omega1 + omega2)/2

dt1 = 0.264
dt2 = 0.287
solar_mass = 5.0289921396853E-31
to_AU = 0.000000000006684587

def omegadot(o1, o2, t1, t2):
    return (-o1+o2)/(t2-t1)

N = 32/5
eta = 1/4
print("---Page 39 ---")
print("Omega = {} rad/s".format(omega))
chirp_mass = c**3/(3**(3/5)*G)*(omegadot(omega1, omega2, dt1, dt2)*omega**(-11/3))**(3/5) * solar_mass*(N*eta)**(-3/5)
print("chirp mass = {} sun mass".format(chirp_mass))


#s_40
w = 2*np.pi*15 #hz
rs_em = 2*G*340/solar_mass/c**2 # 340 solar masses
rs_gr = 2*G*70/solar_mass/c**2 # 70 solar masses
r = ((w**2*2)/(c**2*rs_em))**(-1/3)
r_g = ((w**2*2)/(c**2*rs_gr))**(-1/3)
print("---Page 40 ---")
print("EM-model predicts {:.3e} m, {:.4e} sun diameter from another".format(r, r/1.392e9))
print("GR-model predicts {:.3e} m, {:.4e} sun diameter from another".format(r_g, r_g/1.392e9))

#s_42
Distance = 1402472423.9089 #Ly, 430 Mpc
Milky_Way_radii = 100000#Ly

"""
THe collision happened about 1.4 billion years ago, at a distance of 7000 milky ways
away from us. The distance is larger than superclusters, of which the milky way
in the local group is a member of.
"""

#s_43
omega_start = 15*2*np.pi
lil_r = 1e6
R = G*eta*70/solar_mass*omega_start**2*lil_r**2/( 2*c**3*10 )
print("---Page 43 ---")
print("{} AU".format(R*to_AU))
print("Schwartz radii {:.3e} Our distance {:.3f}".format(rs_gr*to_AU,R*to_AU))


V_theta = np.sqrt( 2*G*70/solar_mass/R/to_AU )
print("{} c".format(V_theta/c)) #Funker null
