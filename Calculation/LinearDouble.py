import numpy as np
import sympy as sym
from sympy import symbols


class LinearDouble:
    def __init__(self):
        H, Htet, Hphi = symbols("H theta_H phi_H", real=True)
        Ms = symbols("Ms", real=True)
        Mtet = symbols("theta_M1 theta_M2", real=True)
        Mphi = symbols("phi_M1 phi_M2", real=True)
        A, Ku2, g = symbols("A K_u2 gamma", real=True)
        J = symbols("J", real=True)
        Hk1, Hk2, H2k = symbols("H_k1 H_k2 H2_k", real=True)
        t = symbols("t", real=True)

        Ev1 = -Ms * H * (sym.sin(Htet) * sym.sin(Mtet[0]) * sym.cos(Hphi - Mphi[0]) +
                         sym.cos(Htet) * sym.cos(Mtet[0])) + \
              A * (sym.sin(Mtet[0]) ** 2) + \
              Ku2 * (sym.sin(Mtet[0]) ** 4)

        Ev2 = -Ms * H * (sym.sin(Htet) * sym.sin(Mtet[1]) * sym.cos(Hphi - Mphi[1]) +
                         sym.cos(Htet) * sym.cos(Mtet[1])) + \
              A * (sym.sin(Mtet[1]) ** 2) + \
              Ku2 * (sym.sin(Mtet[1]) ** 4)

        Es = t * Ev1 + t * Ev2 + \
             J * (sym.sin(Mtet[0]) * sym.sin(Mtet[1]) * sym.cos(Mphi[1] - Mphi[0]) +
                  sym.cos(Mtet[0]) * sym.cos(Mtet[1]))

        self.Ev1_eq = Ev1
        self.Ev2_eq = Ev2
        self.Es_eq = Es

        system_vars = [Mtet[0], Mtet[1], Mphi[0], Mphi[1]]
        system_mat = np.zeros((4, 4), dtype=sym.Add)
        for i in range(2):
            for j in range(4):
                system_mat[i][j] = (g / Ms) * (1 / t) * Es.diff(Mphi[i], system_vars[j])

        for i in range(2, 4):
            for j in range(4):
                system_mat[i][j] = -(g / Ms) * (1 / t) * Es.diff(Mtet[i - 2], system_vars[j])

        w = symbols("w", real=True)
        jw = sym.I * w

        for i in range(4):
            system_mat[i][i] += jw * sym.sin(Mtet[i % 2])

        system_mat = sym.Matrix(system_mat)

        # subs for PMA
        system_mat = system_mat.subs({Hphi: 0, Mphi[0]: 0, Mphi[1]: 0})
        system_mat = system_mat.subs({Ku2: Ms * Hk2 / 4})
        system_mat = system_mat.subs({A: Ms * (Hk1 - Hk2) / 2})
        system_mat = system_mat.subs({J: (t / 2) * Ms * (H2k - Hk1)})
        system_mat.simplify()

        self.w_eq = sym.solve(system_mat.det(), w)
        self.f_eq = np.array(self.w_eq)/(2*sym.pi)
        varset = (
            H, Htet,
            Ms, Mtet[0], Mtet[1],
            Hk1, Hk2, H2k,
            t,
            g
        )

        self.f_res_lamb = []
        for eq in self.f_eq:
            self.f_res_lamb.append(sym.lambdify(varset, eq, 'numpy'))

    def freq(
            self,
            H, Htet,
            Ms, Mtet1, Mtet2,
            Hk1, Hk2, H2k,
            t,
            g,
            ind=1
    ):
        return self.f_res_lamb[ind](
            H, Htet,
            Ms, Mtet1, Mtet2,
            Hk1, Hk2, H2k,
            t,
            g,
        )

