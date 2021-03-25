import numpy as np
import sympy
from sympy import symbols, simplify, lambdify
import pandas as pd
from matplotlib import pyplot as plt


mu0 = 1.2566370621219e-6
gn = 2*np.pi*7.622593285e6 #rad/T
gamma=gn*(2)
H_CGS = 3000. #Oe
B_SI = H_CGS*1e-4 #T
H_SI = 1000./(4*np.pi)*H_CGS #A/m

#     rad/(s*T)       (T)      (T)      (T)          (GHz)  (rad^-1)
print(1.76e11*np.sqrt(150e-3*(150e-3 + 1042e3*mu0))*(1e-9/(2*np.pi)))

print(gamma*B_SI*1e-9)

print(28.02)
print(1.76e2/(2*np.pi))

print(gamma*1e-9)
print(927.400968e8 / 1.054571817 * 2e-9 * 0.3)

print(0.3/mu0)
print(3000.*79.5774715459424)