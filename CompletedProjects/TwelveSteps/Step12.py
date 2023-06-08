import numpy
from matplotlib import pyplot, cm
import sympy
from sympy import init_printing
from sympy.utilities.lambdify import lambdify
import math as math
from mpl_toolkits.mplot3d import Axes3D

x_max = 2
nx = 40
x_vec = numpy.linspace(0, x_max, nx)
dx = x_vec[2] - x_vec[1]

y_max = 2
ny = 40
y_vec = numpy.linspace(0, y_max, ny)
dy = y_vec[2] - y_vec[1]

X, Y = numpy.meshgrid(x_vec, y_vec)

dt = .01
nt = 1400

nu = .1
rho = 1
f = 10

nit = 50

u = numpy.zeros([nx, ny])
v = numpy.zeros([nx, ny])
p = numpy.zeros([nx, ny])
b = numpy.zeros([nx, ny])

def build_up_b(rho, dt, dx, dy, u, v, b):

    b[1:-1, 1:-1] = rho * (
                1 / dt * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dx) + (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dy))) - \
                    - ((u[2:, 1:-1] - u[:-2, 1:-1]) / (2 * dx)) ** 2 - \
                    2 * (((v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dy)) * ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dx))) - \
                    ((v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dy)) ** 2

    # periodic on x = 2
    b[-1, 1:-1] = rho*((u[0, 1:-1] - u[-2, 1:-1])/(2*dy) +
                       (v[-1, 2:] + v[-1, 0:-2])/(2*dx))/dt - \
                    ((u[0, 1:-1] - u[-2, 1:-1])/(2*dx))**2 - \
                    2*((u[-1, 2:] - u[-1, 0:-2]/(2*dy))*(v[0, 1:-1] - v[-2, 1:-1]/(2*dx))) - \
                    ((v[-1, 2:] - v[-1, 0:-2])/(2*dy))**2

    # periodic on x = 0
    # 0 = 1, -2 = -1, -1 = 0
    b[0, 1:-1] = rho * ((u[1, 1:-1] - u[-1, 1:-1]) / (2 * dy) +
                         (v[0, 2:] + v[0, 0:-2]) / (2 * dx)) / dt - \
                  ((u[1, 1:-1] - u[-1, 1:-1]) / (2 * dx)) ** 2 - \
                  2 * ((u[0, 2:] - u[0, 0:-2] / (2 * dy)) * (v[1, 1:-1] - v[-1, 1:-1] / (2 * dx))) - \
                  ((v[0, 2:] - v[0, 0:-2]) / (2 * dy)) ** 2

    return b

def pressure_poisson(p, dx, dy, rho, b):
    pn = numpy.empty_like(p)

    for t in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dy ** 2 +
                          (pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dx ** 2) /
                         (2 * (dx ** 2 + dy ** 2)) -
                         dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) * b[1:-1, 1:-1])

        # Periodic BC Pressure @ x = 2
        p[-1, 1:-1] = (((pn[0, 1:-1] + pn[-2, 1:-1]) * dy ** 2 +
                        (pn[-1, 2:] + pn[-1, 0:-2]) * dx ** 2) /
                       (2 * (dx ** 2 + dy ** 2)) -
                       dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) * b[-1, 1:-1])

        # Periodic BC Pressure @ x = 0
        p[0, 1:-1] = (((pn[1, 1:-1] + pn[-1, 1:-1]) * dy ** 2 +
                       (pn[0, 2:] + pn[0, 0:-2]) * dx ** 2) /
                      (2 * (dx ** 2 + dy ** 2)) -
                      dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) * b[0, 1:-1])

        # Wall boundary conditions, pressure
        p[:, -1] = p[:, -2]  # dp/dy = 0 at y = 2
        p[:, 0] = p[:, 1]  # dp/dy = 0 at y = 0

    return p

def channel_flow(u, v, p, b, dx, dy, dt, rho, f, nu, plt):
    for t in range(nt):
        un = u.copy()
        vn = v.copy()
        print(t/nt)

        b = build_up_b(rho, dt, dx, dy, u, v, b)
        p = pressure_poisson(p, dx, dy, rho, b)

        u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                         un[1:-1, 1:-1] * dt / dx *
                         (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                         vn[1:-1, 1:-1] * dt / dy *
                         (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                         dt / (2 * rho * dx) *
                         (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                         nu * (dt / dx ** 2 *
                               (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]) +
                               dt / dy ** 2 *
                               (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2])) +
                         f * dt)

        v[1:-1, 1:-1] = (vn[1:-1, 1:-1] -
                         un[1:-1, 1:-1] * dt / dx *
                         (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                         vn[1:-1, 1:-1] * dt / dy *
                         (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                         dt / (2 * rho * dy) *
                         (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                         nu * (dt / dx ** 2 *
                               (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1]) +
                               dt / dy ** 2 *
                               (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2])))

        # U period x = 2
        u[-1, 1:-1] = (un[-1, 1:-1] - un[-1, 1:-1] * dt / dx *
                       (un[-1, 1:-1] - un[-2, 1:-1]) -
                       vn[-1, 1:-1] * dt / dy *
                       (un[-1, 1:-1] - un[-1, 0:-2]) -
                       dt / (2 * rho * dx) *
                       (p[0, 1:-1] - p[-2, 1:-1]) +
                       nu * (dt / dx ** 2 *
                             (un[0, 1:-1] - 2 * un[-1, 1:-1] + un[-2, 1:-1]) +
                             dt / dy ** 2 *
                             (un[-1, 2:] - 2 * un[-1, 1:-1] + un[-1, 0:-2])) + f * dt)

        # U period x = 0
        # 0 = 1, -2 = -1, -1 = 0 H
        u[0, 1:-1] = (un[0, 1:-1] - un[0, 1:-1] * dt / dx *
                      (un[0, 1:-1] - un[-1, 1:-1]) -
                      vn[0, 1:-1] * dt / dy *
                      (un[0, 1:-1] - un[0, 0:-2]) -
                      dt / (2 * rho * dx) *
                      (p[1, 1:-1] - p[-1, 1:-1]) +
                      nu * (dt / dx ** 2 *
                            (un[1, 1:-1] - 2 * un[0, 1:-1] + un[-1, 1:-1]) +
                            dt / dy ** 2 *
                            (un[0, 2:] - 2 * un[0, 1:-1] + un[0, 0:-2])) + f * dt)

        # V period x = 2
        v[-1, 1:-1] = (vn[-1, 1:-1] - un[-1, 1:-1] * dt / dx *
                       (vn[-1, 1:-1] - vn[-2, 1:-1]) -
                       vn[-1, 1:-1] * dt / dy *
                       (vn[-1, 1:-1] - vn[-1, 0:-2]) -
                       dt / (2 * rho * dy) *
                       (p[-1, 2:] - p[-1, 0:-2]) +
                       nu * (dt / dx ** 2 *
                             (vn[0, 1:-1] - 2 * vn[-1, 1:-1] + vn[-2, 1:-1]) +
                             dt / dy ** 2 *
                             (vn[-1, 2:] - 2 * vn[-1, 1:-1] + vn[-1, 0:-2])))

        # V period x = 0
        # 0 = 1, -2 = -1, -1 = 0
        v[0, 1:-1] = (vn[0, 1:-1] - un[0, 1:-1] * dt / dx *
                      (vn[0, 1:-1] - vn[-1, 1:-1]) -
                      vn[0, 1:-1] * dt / dy *
                      (vn[0, 1:-1] - vn[0, 0:-2]) -
                      dt / (2 * rho * dy) *
                      (p[0, 2:] - p[0, 0:-2]) +
                      nu * (dt / dx ** 2 *
                            (vn[1, 1:-1] - 2 * vn[0, 1:-1] + vn[-1, 1:-1]) +
                            dt / dy ** 2 *
                            (vn[0, 2:] - 2 * vn[0, 1:-1] + vn[0, 0:-2])))

        u[:, 0] = 0
        u[:, -1] = 0
        v[:, 0] = 0
        v[:, -1] = 0


        if(plt == 1):
            fig = pyplot.figure(figsize=(11, 7), dpi=100)
            # plotting the pressure field as a contour
            pyplot.contourf(Y, X, p, alpha=0.5, cmap=cm.viridis)
            pyplot.colorbar()
            # plotting the pressure field outlines
            pyplot.contour(Y, X, p, cmap=cm.viridis)
            # plotting velocity field
            pyplot.quiver(Y[::2, ::2], X[::2, ::2], u[::2, ::2], v[::2, ::2])
            pyplot.xlabel('X')
            pyplot.ylabel('Y')
            pyplot.pause(.001)
            pyplot.close()
            pyplot.show()

    return u, v, p

plt = 1
u, v, p = channel_flow(u, v, p, b, dx, dy, dt, rho, f, nu, plt)

if(plt == 0):
    fig = pyplot.figure(figsize=(11, 7), dpi=100)
    # plotting the pressure field as a contour
    pyplot.contourf(Y, X, p, alpha=0.5, cmap=cm.viridis)
    pyplot.colorbar()
    # plotting the pressure field outlines
    pyplot.contour(Y, X, p, cmap=cm.viridis)
    # plotting velocity field
    pyplot.quiver(Y[::1, ::1], X[::1, ::1], u[::1, ::1], v[::1, ::1])
    pyplot.xlabel('X')
    pyplot.ylabel('Y')
    pyplot.show()







