"""
File : main.py
[MAIN3 - Schrodinger] Feb. - June 2017
A. Khizar, R. Ndiaye, V. Nicol, M. Pecheux
Referent: F. Aviat

Main file that represents the wave function using the TDSE and a explicit euler numerical model.
"""

# useful imports
import numpy as np
import matplotlib.pyplot as plt
from math import exp


# general variables
x_min = -50.
x_max = 50.
delta_x = 0.1
x0 = 0.
sigma = 1
delta_t = 0.5       # time difference between two steps
nb_steps = 100         # nb of steps to compute (so we get the final psi)

# we make the space discrete
x_values = np.arange(x_min, x_max, delta_x)

# number of subdivisions
n = len(x_values)

# declaration of main matrices
laplacian =  -1 / (2 * delta_x**2) * (-2*np.eye(n) + np.diag(np.ones(n - 1), -1) + np.diag(np.ones(n - 1), 1))
V0 = np.zeros((n, n))                                               # null potential
V1 = np.diag([x ** 2 for x in x_values])                            # harmonic potential
V2 = np.diag([0.5 if abs(x) < 2.0 else 0 for x in x_values])        # rectangular potential
V3 = np.diag([abs(x) for x in x_values])                            # absolute potential
# here you choose which potential you add:
H = laplacian + V0

# graphic representation initialization
plt.figure(1)       # we set the current figure

# declaration of initial wave function
k = 1.
#*np.exp(-1j * k * x)
psi = [exp(-(x - x0) ** 2 / (2*sigma ** 2))*np.exp(-1j * k * x) for x in x_values]
psi2 = psi / np.linalg.norm(psi)

# we get the initial psi values squared (|psi|^2 = presence probability)
y_values = [np.absolute(v) ** 2 for v in psi2]
# we plot the result against the x values
plt.plot(x_values, y_values, label='t = 0')

# computation of the requested steps
for i in range(nb_steps):
    # represent the mid-result
    if i % 20 == 0 and i != 0:
        # we get the final psi values squared (|psi|^2 = presence probability)
        y_values = [np.absolute(v) ** 2 for v in psi]
        # we plot the result against the x values
        plt.plot(x_values, y_values, label='t = %.1f' % (delta_t * i))

    # EXPLICIT EULER METHOD
    # we use the format: 1j to specify that we want a complex number with real part = 0 and imaginary part = 1
    # we then do a dot product to multiply the computed matrix and the vector psi
    #psi = np.dot(np.eye(n) - 1j * delta_t * H, psi)
    # normalization of the vector psi
    #psi = psi / np.linalg.norm(psi)

    # IMPLICIT EULER METHOD
    psi = np.dot(np.linalg.inv(np.eye(n) - 1j * delta_t * H), psi)
    # normalization of the vector psi
    psi = psi / np.linalg.norm(psi)

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
