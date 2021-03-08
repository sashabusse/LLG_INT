import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
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


gamma = 5e8
alpha = 0.018
Hext = np.array(pol2cart(2000.0, np.deg2rad(80.0), 0.0))
# Heff = np.array([0.0, 0.0, -2000.0])
Hkeff = 6100.0


def Heff(Hext, M):
    H_abs, H_tet, H_phi = cart2pol(Hext)
    M_abs, M_tet, M_phi = cart2pol(M)
    r = H_abs * np.cos(H_tet - M_tet)
    tet = H_abs * np.sin(H_tet - M_tet)
    # r = -(H_abs*np.cos(H_tet - M_tet) - 0.5*H_abs*(np.cos(M_tet)**2))
    # tet = -(H_abs*np.sin(H_tet - M_tet) + Hkeff*np.cos(M_tet)*np.sin(M_tet))

    r_ort = np.array([np.sin(M_tet) * np.cos(M_phi), np.sin(M_tet) * np.sin(M_phi), np.cos(M_tet)])
    tet_ort = np.array([np.cos(M_tet) * np.cos(M_phi), np.cos(M_tet) * np.sin(M_phi), -np.sin(M_tet)])
    ret = r * r_ort + tet * tet_ort
    return ret  # pol2cart(r, tet, 0)


def model(t, M):
    Heff_loc = Heff(Hext, M)
    Ms = np.linalg.norm(M)
    first = np.cross(M, Heff_loc)
    second = (alpha / Ms) * np.cross(M, np.cross(M, Heff_loc))
    return -(gamma / (1 + alpha ** 2)) * (first + second)


M0 = np.array([100.0, 0, 150.0])

# t = np.linspace(0,500, 1000) #ps
sol = solve_ivp(model, [0.0, 200.0e-12], M0, method='LSODA')

df = pd.DataFrame()
df['t'] = sol.t
df['x'] = sol.y[0]
df['y'] = sol.y[1]
df['z'] = sol.y[2]
# df = df.set_index('t')


df.plot('t', figsize = (11, 6))
#plt.plot(sol.t, sol.y[0]**2 + sol.y[1]**2 + sol.y[2]**2)
plt.grid()
plt.show()
