import numpy as np
import sympy as sym
from sympy import symbols, lambdify
from sympy.vector import CoordSys3D
from Calculation.geometry_util import sym_gradient_mat


def Energy2Heff(energy, varset, Mx, My, Mz):
    Hx, Hy, Hz = -sym_gradient_mat(energy, Mx, My, Mz, simplify=True)

    varset = tuple(varset)
    Hx_l = lambdify(varset, Hx, 'numpy')
    Hy_l = lambdify(varset, Hy, 'numpy')
    Hz_l = lambdify(varset, Hz, 'numpy')

    return (Hx, Hy, Hz), (Hx_l, Hy_l, Hz_l)


class DemagnetizationEnergy:
    def __init__(self):
        S = CoordSys3D('S')
        Mx, My, Mz = symbols("M_x M_y M_z")

        self.E_eq = 2 * sym.pi * (Mz ** 2)

        self.varset = (Mz,)
        self.H_eq, self.H_l = Energy2Heff(self.E_eq, self.varset, Mx, My, Mz)

    def Heff(self, M):
        return np.array([
            self.H_l[0](M[2]),
            self.H_l[1](M[2]),
            self.H_l[2](M[2])
        ])


class StoenerWolfarthEnergy:
    def __init__(self):
        S = CoordSys3D('S')

        Ku_1ord = symbols("K_u^1ord")
        Mx, My, Mz = symbols('M_x M_y M_z')
        M = Mx * S.i + My * S.j + Mz * S.k

        sin_tet_sqr = (Mx ** 2 + My ** 2) / (Mx ** 2 + My ** 2 + Mz ** 2)
        self.E_eq = Ku_1ord * sin_tet_sqr

        self.varset = (Ku_1ord, Mx, My, Mz)
        self.H_eq, self.H_l = Energy2Heff(self.E_eq, self.varset, Mx, My, Mz)

    def Heff(self, M, Ku_1ord):
        return np.array([
            self.H_l[0](Ku_1ord, M[0], M[1], M[2]),
            self.H_l[1](Ku_1ord, M[0], M[1], M[2]),
            self.H_l[2](Ku_1ord, M[0], M[1], M[2])
        ])


class StoenerWolfarth_2nd_order_Energy:
    def __init__(self):
        S = CoordSys3D('S')

        Ku_2ord = symbols("K_u^2ord")
        Mx, My, Mz = symbols('M_x M_y M_z')
        M = Mx * S.i + My * S.j + Mz * S.k

        sin_tet_sqr = (Mx ** 2 + My ** 2) / (Mx ** 2 + My ** 2 + Mz ** 2)
        self.E_eq = Ku_2ord * (sin_tet_sqr ** 2)

        self.varset = (Ku_2ord, Mx, My, Mz)
        self.H_eq, self.H_l = Energy2Heff(self.E_eq, self.varset, Mx, My, Mz)

    def Heff(self, M, Ku_2ord):
        return np.array([
            self.H_l[0](Ku_2ord, M[0], M[1], M[2]),
            self.H_l[1](Ku_2ord, M[0], M[1], M[2]),
            self.H_l[2](Ku_2ord, M[0], M[1], M[2])
        ])


class IEC_Energy:
    def __init__(self):
        S = CoordSys3D('S')

        J, t = symbols("J t")

        Mx1, My1, Mz1 = symbols('M_x1 M_y1 M_z1')
        M1 = Mx1 * S.i + My1 * S.j + Mz1 * S.k

        Mx2, My2, Mz2 = symbols('M_x2 M_y2 M_z2')
        M2 = Mx2 * S.i + My2 * S.j + Mz2 * S.k

        self.E_eq = (J / t) * M1.dot(M2) / (M1.magnitude() * M2.magnitude())

        self.varset = (
            J, t,
            Mx1, My1, Mz1,
            Mx2, My2, Mz2
        )

        self.H_eq, self.H_l = Energy2Heff(self.E_eq, self.varset, Mx1, My1, Mz1)

    def Heff(self, M_for, M_other, J, t_for):
        return np.array([
            self.H_l[0](J, t_for, M_for[0], M_for[1], M_for[2], M_other[0], M_other[1], M_other[2]),
            self.H_l[1](J, t_for, M_for[0], M_for[1], M_for[2], M_other[0], M_other[1], M_other[2]),
            self.H_l[2](J, t_for, M_for[0], M_for[1], M_for[2], M_other[0], M_other[1], M_other[2])
        ])