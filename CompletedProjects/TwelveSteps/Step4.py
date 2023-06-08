import numpy
from matplotlib import pyplot
import sympy
from sympy import init_printing
from sympy.utilities.lambdify import lambdify
import math as math

x, nu, t = sympy.symbols('x nu t')
phi = (sympy.exp(-(x - 4 * t)**2 / (4 * nu * (t + 1))) +
       sympy.exp(-(x - 4 * t - 2 * sympy.pi)**2 / (4 * nu * (t + 1))))

phi_prime = phi.diff(x)

u = -2 * nu * (phi_prime/phi) + 4
u_func = lambdify((t, x, nu), u)

# ------ #
# u_mat[t+1][x] = u_mat[t][x] - (u_mat[t][x] * (dt/dx) * (u_mat[t][x]-u_mat[t][x-1])) + (nu * (dt/(dx*dx)) * (u_mat[t][x+1]-2*u_mat[t][x]+u_mat[t][x-1]))

nu = .07

total_x = 2*numpy.pi
nx = 101
x_vec = numpy.linspace(0, total_x, nx)
dx = x_vec[2]-x_vec[1]

nt = 300
dt = dx*nu
t_vec = numpy.linspace(0, dt*nu, nt)

u_mat = numpy.zeros([len(t_vec), len(x_vec)])

for t in range(0, len(t_vec)):
    for x in range(0, len(x_vec)):
        u_mat[t][x] = u_func(t_vec[t], x_vec[x], nu)

for t in range(0, len(t_vec)-1):
    for x in range(1, len(x_vec)-1):
        u_mat[t + 1][x] = u_mat[t][x] - (u_mat[t][x] * (dt / dx) * (u_mat[t][x] - u_mat[t][x - 1])) + (
                    nu * (dt / (dx * dx)) * (u_mat[t][x + 1] - 2 * u_mat[t][x] + u_mat[t][x - 1]))
    u_mat[t+1][0] = u_mat[t][0] - (u_mat[t][0] * (dt / dx) * (u_mat[t][0] - u_mat[t][-2])) + (
            nu * (dt / (dx * dx)) * (u_mat[t][1] - 2 * u_mat[t][0] + u_mat[t][-2]))
    u_mat[t][-1] = u_mat[t][0]

for t in range(0, len(t_vec)):
    pyplot.figure(figsize=(11, 7), dpi=100)
    pyplot.plot(x_vec, u_mat[t], marker='o', lw=2, label='bruh')
    pyplot.xlim([0, 2 * numpy.pi])
    pyplot.ylim([0, 10])
    pyplot.pause(.001)
    pyplot.close()

pyplot.show()







