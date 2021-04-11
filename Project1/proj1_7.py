import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
#from progress.bar import Bar
import scipy.integrate as integrate

from proj1_class import CosmoProject1



if "__main__" == __name__:
    Omega = CosmoProject1(N=40)
    Omega.read_file("sndata.txt")
    Omega.find_optimal_parameter_w_omega_w_combination()
