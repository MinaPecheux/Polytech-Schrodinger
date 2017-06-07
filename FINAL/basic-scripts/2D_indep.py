"""
File : 2D_indep.py
[MAIN3 - Schrodinger] Feb. - June 2017
A. Khizar, R. Ndiaye, V. Nicol, M. Pecheux
Referent: F. Aviat
Main file that represents the wave functions in 2D using the TISE.
"""
# useful imports
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

# GENERAL VARIABLES DECLARATION --------------------------------------
x_min = -5.
x_max = 5.
lx = x_max - x_min
y_min = -5.
y_max = 5.
ly = y_max - y_min
n = 50     # number of subdivisions of the space
delta_x = lx / n
delta_y = ly / n

# we make the space discrete
x_values = np.arange(x_min, x_max, delta_x)
y_values = np.arange(y_min, y_max, delta_y)
X, Y = np.meshgrid(x_values, y_values)

# list that corresponds to the superior and inferior diagonals
Diag_SupInf = [0 if i % n == 0 else 1 for i in range(1,n*n)]

# declaration of main matrices
laplacian_x =  1 / delta_x**2 * (-2*np.eye(n*n) + np.diag(Diag_SupInf, -1) + np.diag(Diag_SupInf, 1))
laplacian_y =  1 / delta_y**2 * (-2*np.eye(n*n) + np.diag(np.ones(n*n - n), n) + np.diag(np.ones(n*n - n), -n))
laplacian = -0.5 * (laplacian_x + laplacian_y)

V0 = np.zeros((n*n, n*n))                                    # null potential
#V1 = np.asarray([x ** 2 for x in X])                         # harmonic potential
#V2 = np.diag([0.5 if abs(x) < 2.0 else 0 for x in X])        # rectangular potential
# here you choose which potential you add:
H = laplacian + V0
# ----------------------------------------------------------------

# COMPUTATION + GRAPHIC REPRESENTATION ---------------------------------
# we get the eigen vectors of the H matrix
# to get the eigen functions of the H operator
eigen_values, eigen_vectors = np.linalg.eigh(H)

# we create the figures
fig1, fig2 = plt.figure(1), plt.figure(2)
# we set them for 3d representation
ax1, ax2 = fig1.gca(projection='3d'), fig2.gca(projection='3d')

# orbitals that we want to compute
energy = 1
nx, ny = 2, 2

# analytic solution representation: psi^2 = presence probability
Z = (2/sqrt(lx*ly) * np.multiply(np.sin(nx * np.pi * np.asarray([x + lx / 2 for x in X])/lx), np.sin(ny * np.pi * np.asarray([y + ly / 2 for y in Y])/ly))) ** 2

# we plot the result against the x values
ax1.plot_surface(X, Y, Z, linewidth=0, antialiased=True)

# experimental solution representation: psi^2 = presence probability
psi = list(eigen_vectors[:, energy - 1])
Z = []
for i in range(n):
    Z.append([v ** 2 for v in psi[i*n:(i+1)*n]])

ax2.plot_surface(X, Y, Z, linewidth=0, antialiased=True)
# ----------------------------------------------------------------
plt.show()