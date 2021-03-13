import matplotlib.pyplot as plt
import numpy as np

def U(x, om, ow):
    return -(om/x + ow*x**(-(1+3*w)))

def E(x, om, ow):
    return 1 - om - ow + x - x

x = np.linspace(0.001,1,100)
omega_m = 0
omega_w = 3
omega_m2 = 2
omega_w2 = 1
w = -1



plt.plot(x,U(x, omega_m, omega_w), label = 'U')
plt.plot(x,U(x, omega_m2, omega_w2), label = 'U2')
plt.plot(x, E(x , omega_m, omega_w), label = 'E')
plt.plot(x, E(x, omega_m2, omega_w2), label = 'E2')

plt.xlabel("x")
plt.ylabel("E and U [energy]")
plt.legend()
plt.show()
