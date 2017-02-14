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
x_min = -50
x_max = 50
n = 100       # number of subdivisions of the space
delta_x = (x_max - x_min) / n

# we make the space discrete
x_values = np.arange(x_min, x_max, delta_x)

# declaration of main matrices
laplacian = 1/delta_x**2 * -2*np.eye(n) + np.diag(np.ones(n - 1), -1) + np.diag(np.ones(n - 1), 1)
V0 = np.zeros(n)                                # null potential
V1 = np.asarray([x ** 2 for x in x_values])     # harmonic potential
# here you choose which potential you add:
H = laplacian + V1
# we get the eigen vectors of the H matrix
# to get the eigen functions of the H operator
eigen_values, eigen_vectors = np.linalg.eigh(H)

# graphic representation of the first modes
for i in range(0, 4):
    # we get the right psi values = the corresponding eigen vector
    # of the H matrix

    # THE PROBLEM WAS: the eigen vectors are COLUMN VECTORS!
    y_values = eigen_vectors[:, i]           # using the H eigen vectors
    # --------------------------------
    # The following part is necessary for the null potential to get the sinus, because
    # we have an oscillating curve and we need to switch one point out of two to get back the sinus...
    # (I don't know why?)
    # BUT it must removed for the harmonic potential
    inverse_index = next((k for k, v in enumerate(y_values) if v < 0), -1)
    for j in range(inverse_index, len(y_values), 2):
        y_values[j] = -y_values[j]
    # --------------------------------
    # we plot the result against the x values
    plt.figure(1)
    plt.plot(x_values, y_values, label='mode %d' % (i+1))

# we finish up the graphs and show them
plt.legend()
plt.xlabel('Horizontal position')
plt.ylabel('Wave functions')

plt.show()
