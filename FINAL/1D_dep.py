"""
File : 1D_dep.py
[MAIN3 - Schrodinger] Feb. - June 2017
A. Khizar, R. Ndiaye, V. Nicol, M. Pecheux
Referent: F. Aviat

Main file that represents the wave function using the TDSE and numerical methods:
- explicit euler
- implicit euler
- cranck nicholson

User can choose between a demo of the methods, or of the reflection or the tunnel
effects.
"""
# useful imports
import numpy as np
import matplotlib.pyplot as plt
from math import exp

# NUMERIC METHODS FUNCTIONS ------------------------------------------
def euler_exp(y, H, n, delta_t):
    # we use the format: 1j to specify that we want a complex number with real part = 0 and imaginary part = 1
    # we then do a dot product to multiply the computed matrix and the vector
    y = np.dot(np.eye(n) - 1j * delta_t * H, y)
    # normalization of the vector
    y /= np.linalg.norm(y)
    return y

def euler_imp(y, H, n, delta_t):
    # calculation of new vector
    y = np.dot(np.linalg.inv(np.eye(n) + 1j * delta_t * H), y)
    # normalization of the vector
    y /= np.linalg.norm(y)
    return y

def cranck_nicholson(y, H, n, delta_t):
    # calculation of new vector
    y = np.dot(np.linalg.inv(np.eye(n) + 1j * delta_t/2 * H), np.dot(np.eye(n) - 1j * delta_t / 2 * H, y))
    # normalization of the vector
    y /= np.linalg.norm(y)
    return y
# ---------------------------------------------------------------

def showcase_methods():
    # GENERAL VARIABLES DECLARATION --------------------------------------
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
    V = np.zeros((n, n))                                               # null potential
    # here you choose which potential you add:
    H = laplacian + V
    # ----------------------------------------------------------------

    # INITIALIZATION -----------------------------------------------------
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
    psi_func_methods = ['EULER EXP', 'EULER IMP', 'CRANCK-NICHOLSON']

    for f in range(len(psi_func)):
        plt.figure(f)
        plt.plot(x_values, y_values, label='t = 0')
    # ----------------------------------------------------------------

    # COMPUTATION + GRAPHIC REPRESENTATION ---------------------------------
    for i in range(nb_steps):
        # computation of new values
        psi_func[0] = euler_exp(psi_func[0], H, n, delta_t)
        psi_func[1] = euler_imp(psi_func[1], H, n, delta_t)
        psi_func[2] = cranck_nicholson(psi_func[2], H, n, delta_t)
        # representation of several mid-results
        if i % 10 == 0 and i != 0:
            for f in range(len(psi_func)):
                plt.figure(f)
                y_values = [np.absolute(v) ** 2 for v in psi_func[f]]
                plt.plot(x_values, y_values, label='t = %.1f' % (delta_t * i))
    # ----------------------------------------------------------------

    # GRAPHICS FINISH UP -------------------------------------------------
    for f in range(len(psi_func)):
        plt.figure(f)
        plt.title('Wave function %s with delta_t=%.1f, x0=%.1f, sigma=%.1f' % (psi_func_methods[f], delta_t, x0, sigma))
        plt.legend()
        plt.xlabel('Horizontal position')
        plt.ylabel('Wave function')
    # ----------------------------------------------------------------

    # we show the graphs
    plt.show()


def reflection_effect():
    # GENERAL VARIABLES DECLARATION --------------------------------------
    x_min = -30.
    x_max = 30.
    delta_x = 0.2
    x0 = 10.
    sigma = 1.5
    delta_t = 0.4       # time difference between two steps
    nb_steps = 60     # nb of steps to compute (so we get the final psi)

    # we make the space discrete
    x_values = np.arange(x_min, x_max, delta_x)

    # number of subdivisions
    n = len(x_values)

    # declaration of main matrices
    laplacian = -1 / (2 * delta_x**2) * (-2*np.eye(n) + np.diag(np.ones(n - 1), -1) + np.diag(np.ones(n - 1), 1))
    V = np.zeros((n, n))                                               # null potential
    # here you choose which potential you add:
    H = laplacian + V
    # ----------------------------------------------------------------

    # INITIALIZATION -----------------------------------------------------
    # declaration of initial wave function
    k = -2.
    psi = [exp(-(x - x0) ** 2 / (2*sigma ** 2))*np.exp(-1j * k * x) for x in x_values]

    # we get the initial psi NORMALIZED values squared (|psi|^2 = presence probability)
    psi2 = psi / np.linalg.norm(psi)
    y_values = [np.absolute(v) ** 2 for v in psi2]

    # COMPUTATION + GRAPHIC REPRESENTATION ---------------------------------
    # use EULER IMPLICIT method to compute the wave function movement
    # base gaussian representation
    plt.plot(x_values, y_values, label='t = 0')
    # steps computation and representation
    for i in range(nb_steps):
        psi = euler_imp(psi, H, n, delta_t)
        if i % 15 == 0 and i != 0:
            y_values = [np.absolute(v) ** 2 for v in psi]
            plt.plot(x_values, y_values, label='t = %.1f' % (delta_t * i))
    # ----------------------------------------------------------------

    # GRAPHIC FINISH UP -------------------------------------------------
    plt.title('Reflection effect (with Euler Implicit method)')
    plt.legend()
    plt.xlabel('Horizontal position')
    plt.ylabel('Wave function')
    # ----------------------------------------------------------------

    # we show the graphs
    plt.show()


def tunnel_effect():
    # GENERAL VARIABLES DECLARATION --------------------------------------
    x_min = -50.
    x_max = 50.
    delta_x = 0.2
    x0 = -5.
    sigma = 1.5
    delta_t = 0.4       # time difference between two steps
    nb_steps = 120     # nb of steps to compute (so we get the final psi)

    # we make the space discrete
    x_values = np.arange(x_min, x_max, delta_x)

    # number of subdivisions
    n = len(x_values)

    # declaration of main matrices
    laplacian = -1 / (2 * delta_x**2) * (-2*np.eye(n) + np.diag(np.ones(n - 1), -1) + np.diag(np.ones(n - 1), 1))
    potential_value = 0.5
    V = np.diag([potential_value if abs(x_values[i]) < 1.0 else 0.0 for i in range(n)])                             # potential barrier
    # here you choose which potential you add:
    H = laplacian + V
    # ----------------------------------------------------------------

    # INITIALIZATION -----------------------------------------------------
    # declaration of initial wave function
    k = -1.
    psi = [exp(-(x - x0) ** 2 / (2*sigma ** 2))*np.exp(-1j * k * x) for x in x_values]

    # we get the initial psi NORMALIZED values squared (|psi|^2 = presence probability)
    psi2 = psi / np.linalg.norm(psi)
    y_values = [np.absolute(v) ** 2 for v in psi2]

    # show the potential
    plt.axes(xlim=(x_min, x_max), ylim=(0, 0.1))
    plt.plot(x_values, np.diagonal(V), label='Potential')

    # COMPUTATION + GRAPHIC REPRESENTATION ---------------------------------
    # use EULER IMPLICIT method to compute the wave function movement
    # base gaussian representation
    plt.plot(x_values, y_values, label='t = 0')
    # steps computation and representation
    for i in range(nb_steps):
        psi = euler_imp(psi, H, n, delta_t)
        if i % 30 == 0 and i != 0:
            y_values = [np.absolute(v) ** 2 for v in psi]
            plt.plot(x_values, y_values, label='t = %.1f' % (delta_t * i))
    # ----------------------------------------------------------------

    # GRAPHIC FINISH UP -------------------------------------------------
    plt.title('Tunnel effect (with Euler Implicit method), potential = %.1f' % potential_value)
    plt.legend()
    plt.xlabel('Horizontal position')
    plt.ylabel('Wave function')
    # ----------------------------------------------------------------

    # we show the graphs
    plt.show()

tunnel_effect()
