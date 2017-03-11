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
x_min = -8.
x_max = 8.
n = 100              # number of subdivisions of the space
delta_x = (x_max - x_min) / n
x0 = -1.
sigma = 0.5
delta_t = 0.3       # time difference between two steps
nb_steps = 3        # nb of steps to compute (so we get the final psi)

# we make the space discrete
x_values = np.arange(x_min, x_max, delta_x)

# declaration of main matrices
laplacian = - 1 / delta_x**2 * (-2*np.eye(n) + np.diag(np.ones(n - 1), -1) + np.diag(np.ones(n - 1), 1))
V0 = np.zeros((n, n))                                               # null potential
V1 = np.diag([x ** 2 for x in x_values])                            # harmonic potential
V2 = np.diag([0.5 if abs(x) < 2.0 else 0 for x in x_values])        # rectangular potential
V3 = np.diag([abs(x) for x in x_values])                            # absolute potential
# here you choose which potential you add:
H = laplacian + V0

# graphic representation initialization
plt.figure(1)       # we set the current figure

# declaration of initial wave function
psi = [exp(-(x - x0) ** 2 / (2*sigma ** 2)) for x in x_values]
# we get the initial psi values squared (|psi|^2 = presence probability)
y_values = [np.absolute(v) ** 2 for v in psi]
# we plot the result against the x values
plt.plot(x_values, y_values, label='t = 0')

# computation of the requested steps
for i in range(nb_steps):
    # we use the format: 1j to specify that we want a complex number with real part = 0 and imaginary part = 1
    # we then do a dot product to multiply the computed matrix and the vector psi
    psi = np.dot((np.eye(n) - 1j * delta_t * H), psi)
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
