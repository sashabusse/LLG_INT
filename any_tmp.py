import mpmath as mp
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import sympy as sym


m1 = np.array([[sym.symbols("a11"),sym.symbols("a12")],[sym.symbols("a21"),sym.symbols("a22")]])
m2 = np.array([[1, 1],[1, 1]])

m3 = m1*sym.I

print(m3)