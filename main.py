"""
File : main.py
[MAIN3 - Schrodinger] Feb. - June 2017
A. Khizar, R. Ndiaye, V. Nicol, M. Pecheux
Referent: F. Aviat

Main file that represents the wave function using the TISE.
"""

# useful imports
import numpy as np
import matplotlib.pyplot as plt
from math import sin, pi, sqrt


# general variables
x_min = -5.0
x_max = 5.0
n = 100       # number of subdivisions of the space
delta_x = (x_max - x_min) / n
h = 6.64E-34 / (2*pi)       # Dirac's constant
m = 1E-39

# we make the space discrete
x_values = np.arange(x_min, x_max, delta_x)

# declaration of main matrices
laplacian = - 1 / delta_x**2 * (-2*np.eye(n) + np.diag(np.ones(n - 1), -1) + np.diag(np.ones(n - 1), 1))
V0 = np.zeros((n, n))                                               # null potential
V1 = np.diag([x ** 2 for x in x_values])                            # harmonic potential
V2 = np.diag([0.5 if abs(x) < 2.0 else 0 for x in x_values])        # rectangular potential
V3 = np.diag([abs(x) for x in x_values])                            # absolute potential
# here you choose which potential you add:
H = laplacian + V2
# we get the eigen vectors of the H matrix
# to get the eigen functions of the H operator
eigen_values, eigen_vectors = np.linalg.eigh(H)

# energies that we want to compute
energies = range(0, 4)

# graphic representation of the first energies
# (not separated)
plt.figure(1)       # we set the current figure
for i in energies:
    # we get the right psi values = the corresponding eigen vector
    # of the H matrix (reminder: these are column vectors)
    y_values = [v ** 2 for v in eigen_vectors[:, i]]
    # we plot the result against the x values
    plt.plot(x_values, y_values, label='energy level %d' % (i+1))

# we finish up the graph
plt.title('Wave functions given the energy level')
plt.legend()
plt.xlabel('Horizontal position')
plt.ylabel('Wave functions')

# graphic representation of the first energies
# (separated by arbitrary offset)
plt.figure(2)       # we set the current figure
offset = 0.5
for i in energies:
    # we get the right psi values = the corresponding eigen vector
    # of the H matrix (reminder: these are column vectors)
    y_values = [i*offset + v ** 2 for v in eigen_vectors[:, i]]
    # we plot the result against the x values
    plt.plot(x_values, y_values, label='energy level %d' % (i+1))

# we finish up the graph
plt.title('Wave functions given the energy level with an arbitrary offset = %.1f' % offset)
plt.legend()
plt.xlabel('Horizontal position')
plt.ylabel('Wave functions')

# we show the graphs
plt.show()
