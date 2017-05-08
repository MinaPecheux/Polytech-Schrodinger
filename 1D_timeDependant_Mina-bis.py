"""
File : main.py
[MAIN3 - Schrodinger] Feb. - June 2017
A. Khizar, R. Ndiaye, V. Nicol, M. Pecheux
Referent: F. Aviat

Main file that represents the wave function using the TDSE and numerical methods:
- explicit euler
- implicit euler
- cranck nicholson
- RK 2
- RK 4
"""

# useful imports
import numpy as np
import matplotlib.pyplot as plt
from math import exp

# NUMERIC METHODS FUNCTIONS ------------------------------------------------------------------
def euler_exp(y, H):
    # we use the format: 1j to specify that we want a complex number with real part = 0 and imaginary part = 1
    # we then do a dot product to multiply the computed matrix and the vector
    y = np.dot(np.eye(n) - 1j * delta_t * H, y)
    # normalization of the vector
    y /= np.linalg.norm(y)
    return y


def euler_imp(y, H):
    # calculation of new vector
    y = np.dot(np.linalg.inv(np.eye(n) + 1j * delta_t * H), y)
    # normalization of the vector
    y /= np.linalg.norm(y)
    return y


def cranck_nicholson(y, H):
    # calculation of new vector
    y = np.dot(np.linalg.inv(np.eye(n) + 1j * delta_t/2 * H), np.dot(np.eye(n) - 1j * delta_t / 2 * H, y))
    # normalization of the vector
    y /= np.linalg.norm(y)
    return y 


def rk_2(y, t0, delta_t, H):
    # calculation of new vector
    y = y + delta_t*np.dot(-1j*H, y + delta_t/2*np.dot(-1j*H, y))
    # normalization of the vector
    y /= np.linalg.norm(y)
    return y


def rk_4(y, delta_t, H):
    # get temporary variables
    y2 = y + delta_t/2 * np.dot(-1j*H, y)
    y3 = y + delta_t/2 * np.dot(-1j*H, y2)
    y4 = y + delta_t * np.dot(-1j*H, y3)
    # calculation of new vector
    y = y + delta_t*(1/6*np.dot(-1j*H, y) + 2/6*np.dot(-1j*H, y2) + 2/6*np.dot(-1j*H, y3) + 1/6*np.dot(-1j*H, y4))
    # normalization of the vector
    y /= np.linalg.norm(y)
    return y
# --------------------------------------------------------------------------------------------

# GENERAL VARIABLES DECLARATION --------------------------------------------------------------
x_min = -30.
x_max = 30.
delta_x = 0.2
x0 = 2.
sigma = 1
delta_t = 0.1       # time difference between two steps
nb_steps = 40     # nb of steps to compute (so we get the final psi)

# we make the space discrete
x_values = np.arange(x_min, x_max, delta_x)

# number of subdivisions
n = len(x_values)

# declaration of main matrices
laplacian = -1 / (2 * delta_x**2) * (-2*np.eye(n) + np.diag(np.ones(n - 1), -1) + np.diag(np.ones(n - 1), 1))
V0 = np.zeros((n, n))                                               # null potential
V1 = np.diag([x ** 2 for x in x_values])                            # harmonic potential
V2 = np.diag([0.5 if abs(x) < 2.0 else 0 for x in x_values])        # rectangular potential
V3 = np.diag([abs(x) for x in x_values])                            # absolute potential
# here you choose which potential you add:
H = laplacian + V0
# --------------------------------------------------------------------------------------------

# INITIALIZATION -----------------------------------------------------------------------------
# declaration of initial wave function
k = 0.
psi = [exp(-(x - x0) ** 2 / (2*sigma ** 2))*np.exp(-1j * k * x) for x in x_values]

# we get the initial psi NORMALIZED values squared (|psi|^2 = presence probability)
psi2 = psi / np.linalg.norm(psi)
y_values = [np.absolute(v) ** 2 for v in psi2]

# prepare psi functions array (one per method), all initialized to
# initial psi func (gaussian)
psi_func = []
psi_func.append(psi)    # euler explicit
psi_func.append(psi)    # euler implicit
psi_func.append(psi)    # Cranck-Nicholson
psi_func.append(psi)    # middle point
psi_func.append(psi)    # RK4
psi_func_methods = ['EULER EXP', 'EULER IMP', 'CRANCK-NICHOLSON', 'MIDDLE POINT', 'RK4']

for f in range(len(psi_func)):
    plt.figure(f)
    plt.plot(x_values, y_values, label='t = 0')
# --------------------------------------------------------------------------------------------

# COMPUTATION + GRAPHIC REPRESENTATION -------------------------------------------------------
for i in range(nb_steps):
    # computation of new values
    psi_func[0] = euler_exp(psi_func[0], H)
    psi_func[1] = euler_imp(psi_func[1], H)
    psi_func[2] = cranck_nicholson(psi_func[2], H)
    psi_func[3] = rk_2(psi_func[3], i*delta_t, delta_t, H)
    psi_func[4] = rk_4(psi_func[4], delta_t, H)
    # representation of several mid-results
    if i % 10 == 0 and i != 0:
        for f in range(len(psi_func)):
            plt.figure(f)
            y_values = [np.absolute(v) ** 2 for v in psi_func[f]]
            plt.plot(x_values, y_values, label='t = %.1f' % (delta_t * i))
# --------------------------------------------------------------------------------------------

# GRAPHICS FINISH UP -------------------------------------------------------------------------
for f in range(len(psi_func)):
    plt.figure(f)
    plt.title('Wave function %s with delta_t=%.1f, x0=%.1f, sigma=%.1f' % (psi_func_methods[f], delta_t, x0, sigma))
    plt.legend()
    plt.xlabel('Horizontal position')
    plt.ylabel('Wave function')
# --------------------------------------------------------------------------------------------

# we show the graphs
plt.show()
