import numpy
from matplotlib import pyplot, cm
import sympy
from sympy import init_printing
from sympy.utilities.lambdify import lambdify
import math as math
from mpl_toolkits.mplot3d import Axes3D

def plot2D(x, y, p):
    fig = pyplot.figure(figsize=(11, 7), dpi=100)
    ax = fig.gca(projection='3d')
    X, Y = numpy.meshgrid(x, y)
    surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap=cm.viridis,
            linewidth=1, antialiased=False)
    ax.view_init(30, 225)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')

nt = 100
nx = 100
ny = 100
total_x = 2
total_y = 1

x_vec = numpy.linspace(0, total_x, nx)
dx = x_vec[2]-x_vec[1]

y_vec = numpy.linspace(0, total_y, ny)
dy = y_vec[2]-y_vec[1]

p = numpy.zeros([len(x_vec), len(y_vec)])
pb = numpy.zeros([len(x_vec), len(y_vec)])
b = numpy.zeros([len(x_vec), len(y_vec)])

# boundary conditions
for x in range(0, len(x_vec)):
    for y in range(0, len(y_vec)):
        p[0][y] = 0
        p[-1][y] = 0
        p[x][0] = 0
        p[x][-1] = 0

# source
b[int(nx/4)][int(ny/4)] = 100
b[int(3*nx/4)][3*int(ny/4)] = -100
b[int(nx/2)][int(ny/2)] = -40


for t in range(0, nt):
    pb = p.copy()
    for x in range(1, len(x_vec)-1):
        for y in range(1, len(y_vec)-1):
            p[x][y] = ((p[x+1][y]+p[x-1][y]) * dy**2 + (p[x][y+1]+p[x][y-1]) * dx**2 - b[x][y] * dx**2 * dy**2)/(2*(dx**2 + dy**2))

    for x in range(0, len(x_vec)):
        for y in range(0, len(y_vec)):
            p[0][y] = 0
            p[-1][y] = 0
            p[x][0] = 0
            p[x][-1] = 0

plot2D(x_vec, y_vec, p)
pyplot.show()