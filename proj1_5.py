import numpy as np
import matplotlib.pyplot as plt

c = 1

a_0 = 1

def a(t):
    return a_0*np.exp(np.sqrt(lamda/3)*(t-t0))#Something


def r(a, t):
    return c*np.trapz(1/a,t)



def dL(z,a,t):
    return a_0*(1+z)*r(a,t)

def dL_ana(z,t0,te):
    return a_0*(1+z)*np.sqrt(3)*(np.exp(np.sqrt(1/3)*(t0-te)) - 1)

#Example 1
t0 = 10
te = 1
lamda = 1
t = np.linspace(te,t0,100)
z = np.linspace(0.5,4,100)


plt.plot(z,dL(z,a(t),t),label="Num")
plt.plot(z,dL_ana(z,t0,te),label="Ana")
plt.legend()
plt.show()

#Example 2

z = 2

"""
New a(t) and z:
"""
