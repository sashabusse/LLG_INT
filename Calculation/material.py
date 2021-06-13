import numpy as np
import json
import os
import errno
import Calculation.Constants as Constants


def check_directory(filename):
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


def Hk2Ku(Hk1, Hk2, Ms):
    Ku_2ord = Hk2 * Ms / 4.
    Ku_1ord = (Hk1 * Ms / 2. - 2 * Ku_2ord) + 2 * np.pi * (Ms ** 2)
    return Ku_1ord, Ku_2ord


def Hk2KuJ(Hk1, Hk2, HkJ, Ms, t):
    Ku_1ord, Ku_2ord = Hk2Ku(Hk1, Hk2, Ms)
    J = (t/2.)*Ms*(HkJ - Hk1)
    return Ku_1ord, Ku_2ord, J


def Ku2Hk(Ku_1ord, Ku_2ord, Ms):
    Hk2 = 4 * Ku_2ord / Ms
    Hk1 = 2 * Ku_1ord / Ms + 4 * Ku_2ord / Ms - 4 * np.pi * Ms
    return Hk1, Hk2


def KuJ2Hk(Ku_1ord, Ku_2ord, J, Ms, t):
    Hk1, Hk2 = Ku2Hk(Ku_1ord, Ku_2ord, Ms)
    HkJ = Hk1 + 2*J/(t*Ms)
    return Hk1, Hk2, HkJ


class StaticMaterial:
    def __init__(self):
        self.Ms = 610
        self.gamma = 2.11*Constants.g0_CGS
        self.alpha = 0.1
        self.t = 12e-7
        self.name = 'no name'

        #recalculated parameters
        self.Ku_1ord = 0
        self.Ku_2ord = 0

        Hk1 = 7.5e3
        Hk2 = 6.5e3

        self.update_with_Hk(Hk1, Hk2)

    def update_with_Hk(self, Hk1, Hk2):
        self.Ku_1ord, self.Ku_2ord = Hk2Ku(Hk1, Hk2, self.Ms)

    @classmethod
    def from_Hk(cls, Hk1=7.5e3, Hk2=6.5e3, Ms = 610.):
        inst = cls()

        inst.Ms = Ms
        inst.update_with_Hk(Hk1, Hk2)

        return inst

    @classmethod
    def from_Ku(cls, Ku_1ord, Ku_2ord, Ms):
        inst = cls()

        inst.Ms = Ms
        inst.Ku_1ord = Ku_1ord
        inst.Ku_2ord = Ku_2ord

        return inst

    @classmethod
    def from_dict(cls, d):
        inst = cls()

        for key in d.keys():
            if key == 'Ms':
                inst.Ms = d[key]
            if key == 'gamma':
                inst.gamma = d[key]
            if key == 'alpha':
                inst.alpha = d[key]
            if key == 't':
                inst.t = d[key]
            if key == 'name':
                inst.name = d[key]

            if key == 'Ku_1ord':
                inst.Ku_1ord = d[key]
            if key == 'Ku_2ord':
                inst.Ku_2ord = d[key]

        return inst

    def to_dict(self):
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

    def Hk1_Hk2(self):
        Hk1, Hk2 = Ku2Hk(self.Ku_1ord, self.Ku_2ord, self.Ms)
        return Hk1, Hk2

    def Hk1(self):
        Hk1, Hk2 = self.Hk1_Hk2()
        return Hk1

    def Hk2(self):
        Hk1, Hk2 = self.Hk1_Hk2()
        return Hk2


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

    @classmethod
    def from_J(cls, layer1, layer2, J):
        return cls(layer1, layer2, J)

    def update_with_Hk(self, Hk1, Hk2, HkJ, Ms=None, t=None, zeroJ=False):
        if not (t is None):
            self.l1.st.t = t
            self.l2.st.t = t
        if not (Ms is None):
            self.l1.st.Ms = Ms
            self.l2.st.Ms = Ms

        Ku_1ord, Ku_2ord, self.J = Hk2KuJ(Hk1, Hk2, HkJ, self.l1.st.Ms, self.l1.st.t)
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

