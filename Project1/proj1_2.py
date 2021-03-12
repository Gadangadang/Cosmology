import matplotlib.pyplot as plt
import numpy as np

def U(x, om, ow):
    return -(om/x + ow*x**2)

def E(x, om, ow):
    return 1 - om - ow + x - x

x = np.linspace(0.001,1,100)
omega_m = 0.24
omega_w = 0.64
omega_m2 = 1.1
omega_w2 = 1.5

ediff = E(x , omega_m, omega_w) - U(x, omega_m, omega_w)
ediff2 = E(x , omega_m2, omega_w2) - U(x, omega_m2, omega_w2)

if ediff.any() < 0:
    print("This combo is bad for E")
else:
    print("Is all good for E")

if  ediff2.any() < 0:
    print("This combo is bad for E2")
else:
    print("Is all good for E2")


plt.plot(x,U(x, omega_m, omega_w), label = 'U')
plt.plot(x,U(x, omega_m2, omega_w2), label = 'U2')
plt.plot(x, E(x , omega_m, omega_w), label = 'E')
plt.plot(x, E(x, omega_m2, omega_w2), label = 'E2')

plt.xlabel("x")
plt.ylabel("E and U [energy]")
plt.legend()
plt.show()
