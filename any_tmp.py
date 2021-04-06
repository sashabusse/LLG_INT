import sympy as sym
from sympy import symbols, simplify, lambdify
from sympy.vector import CoordSys3D, matrix_to_vector
import numpy as np
import pandas as pd

def pol2cart(r, tet, phi):
    x = r*np.sin(tet)*np.cos(phi)
    y = r*np.sin(tet)*np.sin(phi)
    z = r*np.cos(tet)
    return x,y,z

def cart2pol(vec):
    x = vec[0]
    y = vec[1]
    z = vec[2]
    r = np.linalg.norm(vec)
    tet = np.arccos(z/r)
    phi = np.arctan(y/x)
    return r, tet, phi

def r_ort(tet, phi):
    x = np.sin(tet)*np.cos(phi)
    y = np.sin(tet)*np.sin(phi)
    z = np.cos(tet)
    return np.array([x, y, z])

def tet_ort(tet, phi):
    x = np.cos(tet)*np.cos(phi)
    y = np.cos(tet)*np.sin(phi)
    z = -np.sin(tet)
    return np.array([x, y, z])

def phi_ort(tet, phi):
    x = -np.sin(phi)
    y = np.cos(phi)
    z = 0
    return np.array([x, y, z])

S = CoordSys3D('S')

Ku = symbols("Ku")
Mx, My, Mz = symbols('Mx My Mz')
M = Mx * S.i + My * S.j + Mz * S.k

sin_tet_sqr = (Mx ** 2 + My ** 2) / (Mx ** 2 + My ** 2 + Mz ** 2)

E = Ku * sin_tet_sqr

Hx = simplify(-sym.diff(E, Mx))
Hy = simplify(-sym.diff(E, My))
Hz = simplify(-sym.diff(E, Mz))

lambda_var_set = (Ku, Mx, My, Mz)
Heff_x = lambdify(lambda_var_set, Hx, 'numpy')
Heff_y = lambdify(lambda_var_set, Hy, 'numpy')
Heff_z = lambdify(lambda_var_set, Hz, 'numpy')




Kuv = 6e5
tet = np.deg2rad(67)
Mv = pol2cart(800e3, tet, 0.)
Hmod = Kuv*np.sin(2*tet)/np.linalg.norm(Mv)


print(Hmod * tet_ort(tet, 0))
print([
    Heff_x(Kuv, Mv[0], Mv[1], Mv[2]),
    Heff_y(Kuv, Mv[0], Mv[1], Mv[2]),
    Heff_z(Kuv, Mv[0], Mv[1], Mv[2]),
])

sm = Mv + \
     2

