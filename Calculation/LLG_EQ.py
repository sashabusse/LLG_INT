import numpy as np
import sympy as sym
from sympy import symbols, simplify, lambdify
from sympy.vector import CoordSys3D


class LLG_EQ:
    def __init__(self):
        S = CoordSys3D('S')

        x, y, z = symbols('x y z')
        dMdt = x * S.i + y * S.j + z * S.k

        Mx, My, Mz = symbols('M_x M_y M_z')
        M = Mx * S.i + My * S.j + Mz * S.k

        Hx, Hy, Hz = symbols('H_x H_y H_z')
        H = Hx * S.i + Hy * S.j + Hz * S.k

        g = sym.Symbol('gamma')
        a = sym.Symbol('alpha')
        Ms = sym.Symbol('M_s')

        self.vec_eq = -dMdt - g * M.cross(H) + (a / Ms) * (M.cross(dMdt))
        sol = sym.solve(self.vec_eq.to_matrix(S), dMdt.to_matrix(S))

        sol[x] = simplify(sol[x])
        sol[y] = simplify(sol[y])
        sol[z] = simplify(sol[z])

        self.sol_lamb = dict()
        lambda_var_set = (
            a, g, Ms,
            Mx, My, Mz,
            Hx, Hy, Hz)

        self.sol_lamb['x'] = lambdify(lambda_var_set, sol[x], 'numpy')
        self.sol_lamb['y'] = lambdify(lambda_var_set, sol[y], 'numpy')
        self.sol_lamb['z'] = lambdify(lambda_var_set, sol[z], 'numpy')

    def right_hand(self, Heff, M, layer, t):
        x = self.sol_lamb['x'](
            layer.alpha_l(t), layer.gamma_l(t), np.linalg.norm(M),
            M[0], M[1], M[2], Heff[0], Heff[1], Heff[2]
        )
        y = self.sol_lamb['y'](
            layer.alpha_l(t), layer.gamma_l(t), np.linalg.norm(M),
            M[0], M[1], M[2], Heff[0], Heff[1], Heff[2]
        )
        z = self.sol_lamb['z'](
            layer.alpha_l(t), layer.gamma_l(t), np.linalg.norm(M),
            M[0], M[1], M[2], Heff[0], Heff[1], Heff[2]
        )
        return [x, y, z]


