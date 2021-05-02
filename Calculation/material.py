import numpy as np
import json
import os
import errno


def check_directory(filename):
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


def Hk2KuJ(H1k1, H1k2, H2k, Ms, t):
    Ku_2ord = H1k2*Ms/4.
    Ku_1ord = (H1k1*Ms/2. - 2*Ku_2ord) + 2*np.pi*(Ms**2)
    J = (t/2.)*Ms*(H2k - H1k1)
    return Ku_1ord, Ku_2ord, J


class StaticMaterial:
    def __init__(self, Ms, gamma=5e8, alpha=0.1, Ku_1ord=7e6, Ku_2ord = 0.0, t=12e-7, name='no name'):
        self.Ms = Ms
        self.gamma = gamma
        self.alpha = alpha
        self.Ku_1ord = Ku_1ord
        self.Ku_2ord = Ku_2ord
        self.t = t
        self.name = name

    @classmethod
    def from_dict(cls, d):
        return cls(
            Ms=d['Ms'],
            gamma=d['gamma'], alpha=d['alpha'],
            Ku_1ord=d['Ku_1ord'], Ku_2ord=d['Ku_2ord'],
            t=d['t'],
            name=d['name']
        )

    def to_json(self):
        result = {
            "name": self.name,
            "Ms": self.Ms,
            "Ku_1ord": self.Ku_1ord,
            "Ku_2ord": self.Ku_2ord,
            "t": self.t,
            "alpha": self.alpha,
            "gamma": self.gamma
        }
        return result

    def print(self):
        print("\t"+self.name)
        print("\t\tKu_1ord = {:.2e}".format(self.Ku_1ord))
        print("\t\tKu_2ord = {:.2e}".format(self.Ku_2ord))
        print("\t\tMs = {:.2f}".format(self.Ms))
        print("\t\tt = {:.2e}".format(self.t))
        print("\t\tgamma = {:.2e}".format(self.gamma))
        print("\t\talpha = {:.4f}".format(self.alpha))


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

    def update_with_Hk(self, H1k1, H1k2, H2k, Ms=None, t=None, zeroJ=False):
        if not (t is None):
            self.l1.st.t = t
            self.l2.st.t = t
        if not (Ms is None):
            self.l1.st.Ms = Ms
            self.l2.st.Ms = Ms

        Ku_1ord, Ku_2ord, J = Hk2KuJ(H1k1, H1k2, H2k, self.l1.st.Ms, self.l1.st.t)

        self.J = J
        if zeroJ:
            self.J = 0.

        self.l1.st.Ku_1ord = Ku_1ord
        self.l1.st.Ku_2ord = Ku_2ord

        self.l2.st.Ku_1ord = Ku_1ord
        self.l2.st.Ku_2ord = Ku_2ord

    @classmethod
    def from_file(cls, name):
        filename = "samples/"+name+".json"
        with open(filename, 'r') as file:
            data = json.load(file)
        layer1 = Material(StaticMaterial.from_dict(data['layer_1']))
        layer2 = Material(StaticMaterial.from_dict(data['layer_2']))
        return cls(layer1, layer2, data['J'])

    def save_sample(self, name):
        filename = "samples/"+name+".json"
        check_directory(filename)
        data = {
            "J": self.J,
            "layer_1": self.l1.st.to_json(),
            "layer_2": self.l2.st.to_json()
        }
        with open(filename, 'w') as file:
            json.dump(data, file, indent=4)

    def print(self):
        print("---------------------------------------------------")
        print("sample parametrs:")
        print("\tJ = {:.2f}".format(self.J))
        it = 1
        for layer in (self.l1, self. l2):
            print()
            print("\tlayer{} parameters:".format(it))
            layer.st.print()
            it += 1
        print("---------------------------------------------------")


def print_Hext_parameters(Hext_abs, Hext_tet):
    print("---------------------------------------------------")
    print("External field:")
    print("\tabs(H) = {:.2e}".format(Hext_abs))
    print("\ttet(H) = {:.2f}".format(np.rad2deg(Hext_tet)))
    print("---------------------------------------------------")

