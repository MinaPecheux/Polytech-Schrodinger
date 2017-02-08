import numpy as np
import matplotlib.pyplot as plt
from math import sin, pi

# general variables
x_min = 0
x_max = 100
n = 100       # number of subdivisions

# declaration of main matrices
laplacian = -2*np.eye(n) + np.diag(np.ones(n - 1), -1) + np.diag(np.ones(n - 1), 1)
V = np.zeros(n)     # null potential
H = laplacian + V
# we get the eigen vectors of the H matrix
# to get the eigen functions of the H operator
eigen_values, eigen_vectors = np.linalg.eigh(H)

# we make the space discrete
x_values = np.arange(x_min, x_max, (x_max - x_min) / n)

# graphic representation
for i in range(1, 3):
    # we get the right psi values = the corresponding eigen vector
    # of the H matrix
    """y_values = eigen_vectors[0]
    for j in range(1, len(eigen_vectors)):
        if j < i:
            y_values += eigen_vectors[j]"""

    y_values = eigen_vectors[i]
    #y_values = [sum(x) for x in zip(*eigen_vectors)]
    #y_values = [sin(pi*i*x/x_max) for x in x_values]
    # we plot the psi function against the x values
    plt.plot(x_values, y_values, label='i = %d' % i)

# we finish up the graph and show it
plt.legend()
plt.xlabel('Horizontal position')
plt.ylabel('Wave functions')
plt.show()