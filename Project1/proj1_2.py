import matplotlib.pyplot as plt
import numpy as np

def U(x):
    return -(omega_m/x + omega_w*x**2)

def E(x):
    return 1 - omega_m - omega_w + x - x

x = np.linspace(0.001,1,100)
omega_m = 0.24
omega_w = 0.64

plt.plot(x,U(x), label = 'U')
plt.plot(x, E(x), label = 'E')
plt.title("omega_m0: {}, omega_w0: {}".format(omega_m, omega_w))
plt.xlabel(r"$\Omega_{m0}$")
plt.ylabel(r"$\Omega_{w0}$")
plt.legend()
plt.show()
