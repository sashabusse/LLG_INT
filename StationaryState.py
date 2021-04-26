import numpy as np
from scipy import optimize
import sympy as sym
from sympy import symbols, simplify, lambdify
from sympy.vector import CoordSys3D, matrix_to_vector
from geometry_util import pol2cart

import matplotlib.pyplot as plt


class StationaryState:
    def __init__(self):
        S = CoordSys3D('S')
        Mx1, My1, Mz1 = symbols("M_x1 M_y1 M_z1")
        Mx2, My2, Mz2 = symbols("M_x2 M_y2 M_z2")
        M1 = Mx1 * S.i + My1 * S.j + Mz1 * S.k
        M2 = Mx2 * S.i + My2 * S.j + Mz2 * S.k

        Hx, Hy, Hz = symbols("H_x H_y H_z")
        H = Hx * S.i + Hy * S.j + Hz * S.k

        Ku1_1ord, Ku2_1ord = symbols("K_u1^1ord K_u2^1ord")
        Ku1_2ord, Ku2_2ord = symbols("K_u1^2ord K_u2^2ord")

        J = symbols("J")
        t1, t2 = symbols("t1 t2")

        self.Ze1 = -H.dot(M1)
        self.Ze2 = -H.dot(M2)

        sin_tet_sqr_1 = (Mx1 ** 2 + My1 ** 2) / (Mx1 ** 2 + My1 ** 2 + Mz1 ** 2)
        sin_tet_sqr_2 = (Mx2 ** 2 + My2 ** 2) / (Mx2 ** 2 + My2 ** 2 + Mz2 ** 2)

        self.SW_D1 = Ku1_1ord * sin_tet_sqr_1 + 2 * sym.pi * (Mz1 ** 2)
        self.SW_D2 = Ku2_1ord * sin_tet_sqr_2 + 2 * sym.pi * (Mz2 ** 2)

        self.SW_2ord_1 = Ku1_2ord * (sin_tet_sqr_1 ** 2)
        self.SW_2ord_2 = Ku2_2ord * (sin_tet_sqr_2 ** 2)

        self.Ev1 = self.Ze1 + self.SW_D1 + self.SW_2ord_1
        self.Ev2 = self.Ze2 + self.SW_D2 + self.SW_2ord_2

        self.Es = J * M1.dot(M2) / (M1.magnitude() * M2.magnitude()) + t1 * self.Ev1 + t2 * self.Ev2

        # ------------------------------------------------------------------
        Ze1_varset = (
            Mx1, My1, Mz1,
            Hx, Hy, Hz
        )
        self.Ze1_lamb = lambdify(Ze1_varset, self.Ze1, 'numpy')

        Ze2_varset = (
            Mx2, My2, Mz2,
            Hx, Hy, Hz
        )
        self.Ze2_lamb = lambdify(Ze2_varset, self.Ze2, 'numpy')

        SW_D1_varset = (
            Mx1, My1, Mz1,
            Ku1_1ord
        )
        self.SW_D1_lamb = lambdify(SW_D1_varset, self.SW_D1, 'numpy')

        SW_D2_varset = (
            Mx2, My2, Mz2,
            Ku2_1ord
        )
        self.SW_D2_lamb = lambdify(SW_D2_varset, self.SW_D2, 'numpy')

        SW_2ord_1_varset = (
            Mx1, My1, Mz1,
            Ku1_2ord
        )
        self.SW_2ord_1_lamb = lambdify(SW_2ord_1_varset, self.SW_2ord_1, 'numpy')

        SW_2ord_2_varset = (
            Mx2, My2, Mz2,
            Ku2_2ord
        )
        self.SW_2ord_2_lamb = lambdify(SW_2ord_2_varset, self.SW_2ord_2, 'numpy')
        # ------------------------------------------------------------------

        Ev1_varset = (
            Mx1, My1, Mz1,
            Hx, Hy, Hz,
            Ku1_1ord,
            Ku1_2ord
        )

        self.Ev1_lamb = lambdify(Ev1_varset, self.Ev1, 'numpy')

        Ev2_varset = (
            Mx2, My2, Mz2,
            Hx, Hy, Hz,
            Ku2_1ord,
            Ku2_2ord
        )

        self.Ev2_lamb = lambdify(Ev2_varset, self.Ev2, 'numpy')

        Es_varset = (
            Mx1, My1, Mz1,
            Mx2, My2, Mz2,
            Hx, Hy, Hz,
            Ku1_1ord, Ku2_1ord,
            Ku1_2ord, Ku2_2ord,
            J,
            t1, t2
        )

        self.Es_lamb = lambdify(Es_varset, self.Es, 'numpy')

    def get_angles_optimize(
            self,
            film,
            Hext,
            tol=None,
            x0=np.array([np.deg2rad(1.), 0.0, np.deg2rad(179.), 0.0]),
            bounds=None,
            random=False,
            random_tries=10,
            print_details=False
    ):
        if bounds is None:
            bounds = np.array([(0, np.pi), (0, 2 * np.pi), (0, np.pi), (0, 2 * np.pi)])

        func_l = lambda x: self.funcmin(x, film, Hext)

        best_res = None
        min_E = 1e10

        n=1
        if random:
            n = random_tries

        for i in range(n):
            if random:
                x0 = np.random.uniform(size=(4,)) * np.array([np.pi, 2 * np.pi, np.pi, 2 * np.pi])

            res = optimize.minimize(
                func_l,
                x0=x0,
                bounds=bounds,
                tol=tol
            )

            if func_l(res.x) < min_E:
                best_res = res
                min_E = func_l(res.x)

        if print_details:
            print(best_res)
        return best_res.x

    def funcmin(self, x, film, Hext):
        tet1 = x[0]
        phi1 = x[1]
        M1 = pol2cart(film.l1.st.Ms, tet1, phi1)

        tet2 = x[2]
        phi2 = x[3]
        M2 = pol2cart(film.l2.st.Ms, tet2, phi2)

        return self.Es_lamb(
            M1[0], M1[1], M1[2],
            M2[0], M2[1], M2[2],
            Hext[0], Hext[1], Hext[2],
            film.l1.st.Ku_1ord, film.l2.st.Ku_1ord,
            film.l1.st.Ku_2ord, film.l2.st.Ku_2ord,
            film.J,
            film.l1.st.t, film.l2.st.t
        )

    def get_angles_grid(
            self,
            film,
            Hext,
            first_grid_sz=100,
            grid_sz=10,
            iterations=10
    ):
        ang0_diap = [0., np.pi]
        ang1_diap = [0., np.pi]

        ind = []
        ang0 = []
        ang1 = []

        for it in range(iterations):
            cur_grid_sz = grid_sz
            if it == 0:
                cur_grid_sz = first_grid_sz

            ang0 = np.linspace(ang0_diap[0], ang0_diap[1], cur_grid_sz)
            ang1 = np.linspace(ang1_diap[0], ang1_diap[1], cur_grid_sz)

            Es_val = np.zeros((cur_grid_sz, cur_grid_sz))
            for i in range(len(ang0)):
                for j in range(len(ang1)):
                    M1 = pol2cart(film.l1.st.Ms, ang0[i], 0.)
                    M2 = pol2cart(film.l2.st.Ms, ang1[j], 0.)
                    Es_val[i][j] = self.Es_lamb(
                        M1[0], M1[1], M1[2],
                        M2[0], M2[1], M2[2],
                        Hext[0], Hext[1], Hext[2],
                        film.l1.st.Ku_1ord, film.l2.st.Ku_1ord,
                        film.l1.st.Ku_2ord, film.l2.st.Ku_2ord,
                        film.J,
                        film.l1.st.t, film.l2.st.t
                    )

            ind = np.unravel_index(Es_val.argmin(), Es_val.shape)

            ang0_diap = [ang0[max(ind[0] - 2, 0)], ang0[min(ind[0] + 2, cur_grid_sz - 1)]]
            ang1_diap = [ang1[max(ind[1] - 2, 0)], ang1[min(ind[1] + 2, cur_grid_sz - 1)]]

        res = [ang0[ind[0]], ang1[ind[1]]]
        return np.array(res)

    def plot_energy(
            self,
            film,
            Hext,
            N=1000,

            ang0_diap_grad=np.array([0., 180.]),
            ang1_diap_grad=np.array([0., 180.]),
            sign=1
    ):
        ang0_diap_rad = np.deg2rad(ang0_diap_grad)
        ang1_diap_rad = np.deg2rad(ang1_diap_grad)

        ang1 = np.linspace(ang0_diap_rad[0], ang0_diap_rad[1], N)
        ang2 = np.linspace(ang1_diap_rad[0], ang1_diap_rad[1], N)

        Es_val = np.zeros((N, N))
        for i in range(len(ang1)):
            M1 = pol2cart(film.l1.st.Ms, ang1[i], 0.)
            M2 = pol2cart(film.l2.st.Ms, ang2, 0.)
            Es_val[i] = self.Es_lamb(
                M1[0], M1[1], M1[2],
                M2[0], M2[1], M2[2],
                Hext[0], Hext[1], Hext[2],
                film.l1.st.Ku_1ord, film.l2.st.Ku_1ord,
                film.l1.st.Ku_2ord, film.l2.st.Ku_2ord,
                film.J,
                film.l1.st.t, film.l2.st.t
            )

        fig, ax = plt.subplots()
        fig.set_figheight(6)
        fig.set_figwidth(6)
        ax.imshow(sign * np.log(Es_val - Es_val.min() + 1e-5))
        fig.tight_layout()

        ticks = np.arange(0, N + 1, N // 5)
        ticks[-1] -= 1
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)

        ax.set_xticklabels(np.linspace(ang0_diap_grad[0], ang0_diap_grad[1], 6))
        ax.set_yticklabels(np.linspace(ang1_diap_grad[0], ang1_diap_grad[1], 6))


st_state = StationaryState()