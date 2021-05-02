import numpy as np
import sympy as sym
from sympy.vector import matrix_to_vector


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


def sym_gradient_mat(eq, x_var, y_var, z_var, simplify=False):
    gx = eq.diff(x_var)
    gy = eq.diff(y_var)
    gz = eq.diff(z_var)

    if simplify:
        gx = sym.simplify(gx)
        gy = sym.simplify(gy)
        gz = sym.simplify(gz)

    return np.array([gx, gy, gz])


def sym_gradient_vec(eq, x_var, y_var, z_var, coord_sys, simplify=False):
    gx, gy, gz = sym_gradient_mat(eq, x_var, y_var, z_var, simplify)
    return matrix_to_vector(sym.Matrix([gx, gy, gz]), coord_sys)