"""
File : main.py
[MAIN3 - Schrodinger] Feb. - June 2017
A. Khizar, R. Ndiaye, V. Nicol, M. Pecheux
Referent: F. Aviat

Main file that represents the wave function using the TISE.
"""

# useful imports
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, sin


# general variables
x_min = 0.0
x_max = 10.0
lx = x_max - x_min

y_min = 0.0
y_max = 10.0
ly = y_max - y_min

n = 100     # number of subdivisions of the space
delta_x = lx / n
delta_y = ly / n

h = 6.64E-34

# we make the space discrete
x_values = np.arange(x_min, x_max, delta_x)
y_values = np.arange(y_min, y_max, delta_y)
X, Y = np.meshgrid(x_values, y_values)

# declaration of main matrices
# list that corresponds to the superior and inferior diagonals
Diag_SupInf = [0 if i % n == 0 else 1 for i in range(1,n*n)]

# declaration of main matrices
laplacian_x =  1 / delta_x**2 * (-2*np.eye(n*n) + np.diag(Diag_SupInf, -1) + np.diag(Diag_SupInf, 1))
laplacian_y =  1 / delta_y**2 * (-2*np.eye(n*n) + np.diag(np.ones(n*n - n), n) + np.diag(np.ones(n*n - n), -n))
laplacian = -0.5 * (laplacian_x + laplacian_y)

V0 = np.zeros((n*n, n*n))
V1 = []
for x in x_values:
    for y in y_values:
        V1.append(x**2 + y**2)
V1 = np.diag(V1)                    # harmonic potential
V2 = np.asarray([abs(round(sin(3*x))) for x in x_values])       # rectangular potential
# here you choose which potential you add:
H_x = laplacian_x + V1
H_y = laplacian_y + V1
# we get the eigen vectors of the H matrix
# to get the eigen functions of the H operator
eigen_values_x, eigen_vectors_x = np.linalg.eigh(H_x)
eigen_values_y, eigen_vectors_y = np.linalg.eigh(H_y)

# orbitals that we want to compute
nx, ny = 2,2

# --------------------------------------------------------------
# ANALYTIC solution representation: psi = wave function
# --------------------------------------------------------------
# to be coherent with our experimental solution
# our analytic solution will also be representing the nx+1 and the ny+1 orbitals
Z = 2/sqrt(lx*ly) * np.multiply(np.sin((nx+1) * np.pi * X/lx), np.sin((ny+1) * np.pi * Y/ly))
# we set the figure
plt.figure(1)
# we get the z limits for the intensity ruler
z_min, z_max = np.min(Z), np.max(Z)
# we produce the graph
#plt.pcolor(X, Y, Z, cmap='RdBu', vmin=z_min, vmax=z_max)
plt.pcolor(X, Y, Z, cmap=plt.cm.Oranges, vmin=z_min, vmax=z_max)
plt.colorbar()

# --------------------------------------------------------------
# COMPUTATIONAL solution representation: psi = wave function
# --------------------------------------------------------------
Z = [eigen_vectors_x[:, nx] * y for y in eigen_vectors_y[:, ny]]
# we set the figure
plt.figure(2)
# we get the z limits for the intensity ruler
z_min, z_max = np.min(Z), np.max(Z)
# we produce the graph
plt.pcolor(X, Y, Z, cmap=plt.cm.Oranges, vmin=z_min, vmax=z_max)
plt.colorbar()
plt.title('Wave function in 2-D: intensity\nOrbitals 2,2, null potential')

plt.show()
