import numpy
from matplotlib import pyplot, cm
import sympy
from sympy import init_printing
from sympy.utilities.lambdify import lambdify
import math as math
from mpl_toolkits.mplot3d import Axes3D

sigma = .2
nu = .05

nx = 100
ny = 100
nt = 200

total_x = 5
total_y = 5

x_vec = numpy.linspace(0, total_x, nx)
dx = x_vec[2]-x_vec[1]

y_vec = numpy.linspace(0, total_y, ny)
dy = y_vec[2]-y_vec[1]

dt = sigma * dx * dy / nu
t_vec = numpy.linspace(0, dt*nt, nt)

u_mat = numpy.zeros([len(t_vec), len(x_vec), len(y_vec)])

for t in range(0, len(t_vec)):
    for x in range(0, len(x_vec)):
        for y in range(0, len(y_vec)):
            if (x > nx * 1 / 4 and x < nx * 1 / 2 and y > ny * 1 / 4 and y < ny * 1 / 2):
                u_mat[t][x][y] = 2

for t in range(0, len(t_vec)-1):
    for x in range(1, len(x_vec)-1):
        for y in range(1, len(y_vec)-1):
            u_mat[t+1][x][y] = u_mat[t][x][y] + ((nu*dt)/(dx**2)) * (u_mat[t][x+1][y] - 2*u_mat[t][x][y] + u_mat[t][x-1][y]) + ((nu*dt)/(dy**2)) * (u_mat[t][x][y+1] - 2*u_mat[t][x][y] + u_mat[t][x][y-1])
    print(t/nt)

for t in range(0, len(t_vec)-1):
    '''
    fig = pyplot.figure()
    ax = fig.gca(projection='3d')
    X, Y = numpy.meshgrid(x_vec, y_vec)
    surf = ax.plot_surface(X, Y, u_mat[t], cmap=cm.viridis,
            linewidth=0, antialiased=False)
    pyplot.pause(.0001)
    pyplot.close()
    '''
    pyplot.imshow(u_mat[t])
    pyplot.pause(.0001)
    pyplot.close()

pyplot.show()





