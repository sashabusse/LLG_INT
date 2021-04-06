import numpy as np


class StaticMaterial:
    def __init__(self, Ms, gamma=5e8, alpha=0.1, Ku_1ord=7e6, Ku_2ord = 0.0, t=12e-7):
        self.Ms = Ms
        self.gamma = gamma
        self.alpha = alpha
        self.Ku_1ord = Ku_1ord
        self.Ku_2ord = Ku_2ord
        self.t = t


class Material:
    def __init__(self, static):
        self.st = static

        self.Ms_l = lambda t: self.st.Ms
        self.gamma_l = lambda t: self.st.gamma
        self.alpha_l = lambda t: self.st.alpha
        self.Ku_1ord_l = lambda t: self.st.Ku_1ord
        self.Ku_2ord_l = lambda t: self.st.Ku_2ord
        self.t_l = lambda t: self.st.t

    def set_Ku_1ord_lambda(self, Ku_1ord_lambda):
        self.Ku_1ord_l = Ku_1ord_lambda


class LayeredFilm:
    def __init__(self, layer1, layer2, J):
        self.l1 = layer1
        self.l2 = layer2
        self.J = J


# print parameters
def print_film_parameters(film):
    print("---------------------------------------------------")
    print("film parametrs:")
    print("\tJ = {:.2f}".format(film.J))
    it = 1
    for layer in (film.l1, film.l2):
        print()
        print("\tlayer{} parameters:".format(it))
        print("\t\tKu_1ord{} = {:.2e}".format(it, layer.st.Ku_1ord))
        print("\t\tKu_2ord{} = {:.2e}".format(it, layer.st.Ku_2ord))
        print("\t\tMs{} =".format(it), layer.st.Ms)
        print("\t\tt{} = {:.2e}".format(it, layer.st.t))
        print("\t\tgamma{} = {:.2e}".format(it, layer.st.gamma))
        print("\t\talpha{} =".format(it), layer.st.alpha)
        it += 1

    print("---------------------------------------------------")


def print_Hext_parameters(Hext_abs, Hext_tet):
    print("---------------------------------------------------")
    print("External field:")
    print("\tabs(H) = {:.2e}".format(Hext_abs))
    print("\ttet(H) = {:.2f}".format(np.rad2deg(Hext_tet)))
    print("---------------------------------------------------")
