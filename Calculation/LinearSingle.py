import numpy as np
import sympy as sym
from sympy import symbols


class LinearSingle:
    def __init__(self):
        H, Htet, Hphi = symbols("H theta_H phi_H")
        Ms, Mtet, Mphi = symbols("Ms theta_M phi_M")
        A, Ku1, Ku2, gamma, alpha = symbols("A K_u1 K_u2 gamma alpha")
        Hk1, Hk2 = symbols("H_k1 H_k2")

        E = -Ms * H * (sym.sin(Htet) * sym.sin(Mtet) * sym.cos(Hphi - Mphi) +
                       sym.cos(Htet) * sym.cos(Mtet)) + \
            A * (sym.sin(Mtet) ** 2) + \
            Ku2 * (sym.sin(Mtet) ** 4)

        Etp = E.diff(Mtet).diff(Mphi)
        Ett = E.diff(Mtet).diff(Mtet)
        Epp = E.diff(Mphi).diff(Mphi)

        beta = alpha*gamma/(2*(1-alpha**2)*Ms*(sym.sin(Mtet)**2))*(Epp + Ett*(sym.sin(Mtet)**2))
        w_sqr_K = (gamma / (Ms * sym.sin(Mtet))) ** 2 * (Ett * Epp - Etp ** 2) / (1 + alpha**2)

        beta_sym = symbols("beta")
        w_sqr_full = (gamma / (Ms * sym.sin(Mtet))) ** 2 * (Ett * Epp - Etp ** 2) / (1 + alpha**2) + \
                     (beta*alpha*gamma/Ms)*(Epp+Ett*sym.sin(Mtet)**2)


        w_sqr_K = w_sqr_K.subs({Hphi: 0, Mphi: 0})
        w_sqr_H = w_sqr_K.subs({
            Ku2: Ms * Hk2 / 4,
            A: Ms * (Hk1 - Hk2) / 2
        })
        # w_sqr_H = w_sqr_H.subs({Hk2: 0})

        w_sqr_K = w_sqr_K.subs({
            A: Ku1-2*sym.pi*(Ms**2)
        })

        beta = beta.subs({
            Hphi: 0, Mphi: 0,
            Ku2: Ms * Hk2 / 4,
            A: Ms * (Hk1 - Hk2) / 2
        })
        # tau = tau.subs({Hk2: 0})

        w_sqr_full = w_sqr_full.subs({
            Hphi: 0, Mphi: 0,
            Ku2: Ms * Hk2 / 4,
            A: Ms * (Hk1 - Hk2) / 2
        })
        # w_sqr_full = w_sqr_full.subs({Hk2: 0})

        w_sqr_H = sym.simplify(w_sqr_H)
        w_sqr_K = sym.simplify(w_sqr_K)
        beta = sym.simplify(beta)
        w_sqr_full = sym.simplify(w_sqr_full)

        self.f_H_eq = sym.sqrt(w_sqr_H)/(2*sym.pi)
        self.f_K_eq = sym.sqrt(w_sqr_K)/(2*sym.pi)
        self.beta_eq = beta
        self.f_full_eq = sym.sqrt(w_sqr_full)/(2*sym.pi)

        f_varset = (
            H, Htet,
            Mtet,
            Hk1, Hk2,
            gamma, alpha
        )
        self.f_lamb = sym.lambdify(f_varset, self.f_H_eq, 'numpy')
        self.beta_lamb = sym.lambdify(f_varset, self.beta_eq, 'numpy')

        f_full_varset = (
            H, Htet,
            Mtet,
            Hk1, Hk2,
            gamma, alpha,
            beta_sym
        )
        self.f_full_lamb = sym.lambdify(f_full_varset, self.f_full_eq, 'numpy')

    def freq(
            self,
            H, Htet,
            Mtet,
            Hk1, Hk2,
            gamma, alpha
    ):
        return self.f_lamb(H, Htet, Mtet, Hk1, Hk2, gamma, alpha)

    def tau(
            self,
            H, Htet,
            Mtet,
            Hk1, Hk2,
            gamma, alpha
    ):
        return 1./self.beta_lamb(H, Htet, Mtet, Hk1, Hk2, gamma, alpha)

