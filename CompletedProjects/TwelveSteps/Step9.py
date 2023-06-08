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
            linewidth=0, antialiased=False)
    ax.view_init(30, 225)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')

def laplace2D(p, yr, dx, dy, target):
    norm = 1
    pn = numpy.empty_like(p)
    while norm>target:
        pn = p.copy()
        for x in range(1, len(p)-1):
            for y in range(1, len(p[0])-1):
                p[x][y] = (dy**2*(pn[x+1][y]+pn[x-1][y])+dx**2*(pn[x][y+1]+pn[x][y-1]))/(2*(dx**2+dy**2))

        p[:, 0] = 0  # p = 0 @ x = 0
        p[:, -1] = yr  # p = y @ x = 2
        p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
        p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1

        norm = (numpy.sum(numpy.abs(p[:]) - numpy.abs(pn[:])) /
                numpy.sum(numpy.abs(pn[:])))

    return p

total_x = 2
total_y = 2

nx = 40
ny = 40

x_vec = numpy.linspace(0, total_x, nx)
dx = x_vec[2] - x_vec[1]

y_vec = numpy.linspace(0, total_y, ny)
dy = y_vec[2] - y_vec[1]

p_mat = numpy.zeros([len(x_vec), len(y_vec)])


p_mat[:, 0] = 0  # p = 0 @ x = 0
p_mat[:, -1] = y_vec  # p = y @ x = 2
p_mat[0, :] = p_mat[1, :]  # dp/dy = 0 @ y = 0
p_mat[-1, :] = p_mat[-2, :]  # dp/dy = 0 @ y = 1

plot2D(x_vec, y_vec, laplace2D(p_mat, y_vec, dx, dy, 1e-4))
pyplot.show()
