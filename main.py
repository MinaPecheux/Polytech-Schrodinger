"""
File : main.py
[MAIN3 - Schrödinger] Feb. - June 2017
A. Khizar, R. Ndiaye, V. Nicol, M. Pecheux
Referent: F. Aviat

Main file that represents the wave function using the TISE.
"""

# useful imports
import numpy as np
import matplotlib.pyplot as plt
from math import sin, pi


# general variables
x_min = 0
x_max = 100
n = 100       # number of subdivisions of the space
delta_x = (x_max - x_min) / n

# declaration of main matrices
laplacian = 1/delta_x**2 * -2*np.eye(n) + np.diag(np.ones(n - 1), -1) + np.diag(np.ones(n - 1), 1)
V = np.zeros(n)     # null potential
H = laplacian + V
# we get the eigen vectors of the H matrix
# to get the eigen functions of the H operator
eigen_values, eigen_vectors = np.linalg.eigh(H)

# we make the space discrete
x_values = np.arange(x_min, x_max, delta_x)

# graphic representation of the first modes
for i in range(0, 2):
    # we get the right psi values = the corresponding eigen vector
    # of the H matrix
    y_values = [sin((i+1)*pi*x/x_max) for x in x_values]  # analytic solution
    y_values2 = eigen_vectors[i]                            # using the H eigen vectors
    # we plot the both results against the x values
    plt.figure(1)
    plt.plot(x_values, y_values, label='mode %d' % (i+1))
    plt.figure(2)
    plt.plot(x_values, y_values2, label='mode %d' % (i+1))

# we finish up the graphs and show them
plt.figure(1)
plt.legend()
plt.xlabel('Horizontal position')
plt.ylabel('Wave functions')

plt.figure(2)
plt.legend()
plt.xlabel('Horizontal position')
plt.ylabel('Wave functions')

plt.show()
