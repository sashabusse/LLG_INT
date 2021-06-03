import mpmath as mp
import numpy as np

from itertools import product
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

for i, j in product(a, b):
    print(i)
    print(j)
    print()

