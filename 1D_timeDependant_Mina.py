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


# numeric methods functions
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
    y = np.dot(np.dot(np.linalg.inv(np.eye(n) + 1j * delta_t / 2 * H), np.eye(n) - 1j * delta_t / 2 * H), y)
    # normalization of the vector
    y /= np.linalg.norm(y)
    return y


def rk_2(y, H):
    # calculation of new vector
    A = np.eye(n) - 1j * delta_t / 2 * H
    yNew = y - 1j * delta_t * np.dot(H, np.dot(A, y))
    y = yNew
    # normalization of the vector
    y /= np.linalg.norm(y)
    return y


def rk_4(y, H):
    # calculation of new vector
    A = np.eye(n) - 1j * delta_t / 2 * H
    C = 1j * delta_t * H
    #y = np.dot(np.eye(n) - 1/6*np.dot(C, (4*np.eye(n) + (2+delta_t)*A - C + 1/2*(np.dot(np.dot(C, C), A)))), y)
    E = np.dot(-1j*H, y)
    y = y + delta_t*(1/6*E - np.dot(2*1j*H/6, (y - delta_t / 2 * 1j * np.dot(H, y))) - np.dot(2*1j*H/6, y + delta_t/2*(y - np.dot(1j * H * delta_t/2, y))) - np.dot(1j*H/6, y - delta_t*1j*np.dot(H, y + delta_t/2*np.dot(-1j*H, y + delta_t/2*E))))
    # normalization of the vector
    y /= np.linalg.norm(y)
    return y


# general variables
x_min = -30.
x_max = 30.
delta_x = 0.2
x0 = 2.
sigma = 2
delta_t = 1       # time difference between two steps
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

# graphic representation initialization
plt.figure(1)       # we set the current figure

# declaration of initial wave function
k = 0.
psi = [exp(-(x - x0) ** 2 / (2*sigma ** 2))*np.exp(-1j * k * x) for x in x_values]

# we get the initial psi NORMALIZED values squared (|psi|^2 = presence probability)
psi2 = psi / np.linalg.norm(psi)
y_values = [np.absolute(v) ** 2 for v in psi2]
# we plot the result against the x values
plt.plot(x_values, y_values, label='t = 0')

# computation of the requested steps
for i in range(nb_steps):
    # represent several mid-results
    if i % 10 == 0 and i != 0:
        # we get the psi values squared (|psi|^2 = presence probability)
        y_values = [np.absolute(v) ** 2 for v in psi]
        # we plot the result against the x values
        plt.plot(x_values, y_values, label='t = %.1f' % (delta_t * i))

    # (all methods are declared in functions at the beginning of the file)
    # EXPLICIT EULER METHOD
    #psi = euler_exp(psi, H)

    # IMPLICIT EULER METHOD
    psi = euler_imp(psi, H)

    # CRANCK-NICHOLSON METHOD
    #psi = cranck_nicholson(psi, H)

    # RUNGE-KUTTA 2 = MIDDLE POINT METHOD
    #psi = rk_2(psi, H)

    # RUNGE-KUTTA 4 METHOD
    #psi = rk_4(psi, H)

# we get the final psi values squared (|psi|^2 = presence probability)
y_values = [np.absolute(v) ** 2 for v in psi]
# we plot the result against the x values
plt.plot(x_values, y_values, label='t = %.1f' % (delta_t * nb_steps))

# we finish up the graph
plt.title('Wave function with delta_t=%.1f, x0=%.1f, sigma=%.1f' % (delta_t, x0, sigma))
plt.legend()
plt.xlabel('Horizontal position')
plt.ylabel('Wave function')

# we show the graphs
plt.show()
