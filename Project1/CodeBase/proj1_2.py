import matplotlib.pyplot as plt
import numpy as np

def U(x, om, ow, w):
    return -(om/x + ow*x**(-(1+3*w)))

def E(x, om, ow):
    return 1 - om - ow + x - x

x = np.linspace(0,3,3001)
a0index = int(np.where( np.abs(x - 1.) < 1e-5 )[0])


omega_m = 0.5
omega_w = 3
w = -1

omega_m2 = 2.5
omega_w2 = -0.5
w2 = -1

omega_m3 = 0.3
omega_w3 = 0.64
w3 = -1



#Case 1
plt.plot(x,U(x, omega_m, omega_w, w), label = 'U')
plt.plot(x, E(x, omega_m, omega_w), label = 'E')
plt.plot( 1, U(x, omega_m, omega_w, w)[a0index], "ro", label=r"$a_0$" )
plt.ylim(-30,10)

plt.xlabel(r"x $[a/a_0]$")
plt.ylabel("E and U [energy]")
plt.legend()
#plt.title(r"Case 1: w = {}, $\Omega_{m 0}$ = {}, $\Omega_{w 0}$ = {} ".format(w, omega_m, omega_w))
plt.title("No Big Bang")
plt.savefig("images/case1_2.jpeg")
plt.show()



#Case two
plt.plot(x,U(x, omega_m2, omega_w2, w2), label = 'U')
plt.plot(x, E(x, omega_m2, omega_w2), label = 'E')
plt.plot( 1, U(x, omega_m2, omega_w2, w2)[a0index], "ro", label=r"$a_0$" )
plt.ylim(-30,10)
plt.xlabel(r"x $[a/a_0]$")
plt.ylabel("E and U [energy]")

plt.legend()
#plt.title(r"Case 2: w = {}, $\Omega_{m 0}$ = {}, $\Omega_{w 0}$ = {} ".format(w2, omega_m2, omega_w2))
plt.title("Big Crunch")
plt.savefig("images/case2_2.jpeg")
plt.show()

#Case 3
plt.plot(x,U(x, omega_m3, omega_w3, w3), label = 'U')
plt.plot(x, E(x, omega_m3, omega_w3), label = 'E')
plt.plot( 1, U(x, omega_m3, omega_w3, w3)[a0index], "ro", label=r"$a_0$" )
plt.ylim(-30,10)

plt.xlabel(r"x $[a/a_0]$")
plt.ylabel("E and U [energy]")
#plt.title(r"Case 3: w = {}, $\Omega_{m 0}$ = {}, $\Omega_{w 0}$ = {} ".format(w3, omega_m3, omega_w3))
plt.title("Big Freeze")
plt.legend()
plt.savefig("images/case3_2.jpeg")

plt.show()
