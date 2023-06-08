import numpy
from matplotlib import pyplot, cm
import sympy
from sympy import init_printing
from sympy.utilities.lambdify import lambdify
import math as math
from mpl_toolkits.mplot3d import Axes3D

x_max = 2
y_max = 2
nx = 50
ny = 50
c = 1
dx = x_max / (nx - 1)
dy = y_max / (ny - 1)
x_vec = numpy.linspace(0, x_max, nx)
y_vec = numpy.linspace(0, y_max, ny)
X, Y = numpy.meshgrid(x_vec, y_vec)

rho = 1
nu = .1
dt = .001

def build_up_b(b, dt, dx, dy, u, v, rho):
    b[1:-1, 1:-1] = rho*(1/dt*((u[2:, 1:-1]-u[0:-2, 1:-1])/(2*dx) + (v[1:-1, 2:]-v[1:-1, 0:-2])/(2*dy))) - \
                    - ((u[2:, 1:-1]-u[:-2, 1:-1])/(2*dx))**2 - \
                    2 * (((v[1:-1, 2:] - v[1:-1, 0:-2])/(2*dy)) * ((v[2:, 1:-1] - v[0:-2, 1:-1])/(2*dx))) - \
                    ((v[1:-1, 2:] - v[1:-1, 0:-2])/(2*dy))**2

    return b

def pressure_poisson(p, dx, dy, b, rho):
    pn = numpy.empty_like(p)
    pn = p.copy()

    for t in range(nit):
        pn = p.copy()

        p[1:-1, 1:-1] = ((pn[2:, 1:-1]+pn[0:-2, 1:-1])*(dy**2) + (pn[1:-1, 2:]+pn[1:-1, 0:-2])*(dx**2))/(2*(dx**2 + dy**2)) - \
                        ((rho*(dx**2)*(dy**2))/(2*(dx**2+dy**2))) * b[1:-1, 1:-1]

        p[-1, :] = p[-2, :]  # dp/dx = 0 at x = 2
        p[:, 0] = p[:, 1]  # dp/dy = 0 at y = 0
        p[0, :] = p[1, :]  # dp/dx = 0 at x = 0
        p[:, -1] = 0  # p = 0 at y = 2

       # p[10:20, 10:20] = 0
       # p[30:40, 30:40] = 0

    return p

def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, plt):
    un = numpy.empty_like(u)
    vn = numpy.empty_like(v)
    b = numpy.zeros([nx, ny])

    for t in range(0, nt):
        print(t/nt)
        un = u.copy()
        vn = v.copy()

        b = build_up_b(b, dt, dx, dy, u, v, rho)
        p = pressure_poisson(p, dx, dy, b, rho)

        u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                         un[1:-1, 1:-1] * dt / dx *
                         (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                         vn[1:-1, 1:-1] * dt / dy *
                         (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                         dt / (2 * rho * dx) * (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                         nu * (dt / dx ** 2 *
                               (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]) +
                               dt / dy ** 2 *
                               (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2])))

        v[1:-1, 1:-1] = (vn[1:-1, 1:-1] -
                         un[1:-1, 1:-1] * dt / dx *
                         (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                         vn[1:-1, 1:-1] * dt / dy *
                         (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                         dt / (2 * rho * dy) * (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                         nu * (dt / dx ** 2 *
                               (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1]) +
                               dt / dy ** 2 *
                               (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2])))

        u[:, 0] = 0
        u[0, :] = 0
        u[-1, :] = 0
        u[:, -1] = 1  # set velocity on cavity lid equal to 1

        v[:, 0] = 0
        v[:, -1] = 0
        v[0, :] = 0
        v[-1, :] = 0

        #wall in the middle
        u[10:20, 20:30] = 0
        v[10:20, 20:30] = 0

        u[30:40, 30:40] = 0
        v[30:40, 30:40] = 0

        u[10:15, 10:15] = 0
        v[10:15, 10:15] = 0


        if(plt == 1):
            fig = pyplot.figure(figsize=(11, 7), dpi=100)
            # plotting the pressure field as a contour
            pyplot.contourf(Y, X, p, alpha=0.5, cmap=cm.viridis)
            clb = pyplot.colorbar()
            clb.set_label('Pressure')
            # plotting the pressure field outlines
            pyplot.contour(Y, X, p, cmap=cm.viridis)
            # plotting velocity field
            pyplot.quiver(Y[::2, ::2], X[::2, ::2], u[::2, ::2], v[::2, ::2])
            pyplot.xlabel('X')
            pyplot.ylabel('Y')

            rect1 = pyplot.Rectangle((29 * dx, 29 * dy), 10 * dx, 10 * dy, color='purple')
            rect2 = pyplot.Rectangle((9 * dx, 19 * dy), 10 * dx, 10 * dy, color='purple')
            rect3 = pyplot.Rectangle((29 * dx, 9 * dy), 5 * dx, 5 * dy, color='purple')
            pyplot.gca().add_patch(rect1)
            pyplot.gca().add_patch(rect2)
            pyplot.gca().add_patch(rect3)

            pyplot.pause(.0001)
            pyplot.close()

        if(plt == 2):
            a = numpy.zeros((nx, ny))
            a[1:-1] = u[1:-1]*u[1:-1] + v[1:-1]*v[1:-1]
            pyplot.imshow(a)
            pyplot.pause(.0001)
            pyplot.close()

    return u, v, p

u = numpy.zeros([nx, ny])
v = numpy.zeros([nx, ny])
p = numpy.zeros([nx, ny])

nt = 500
nit = 100

# plot style--
plt = 3 # 1 = time quiver, 2 = time magnitude, 3 = final quiver, 4 = final magnitude, 5 = final pressure im
# ------------

u, v, p = cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, plt)

if(plt == 3):
    fig = pyplot.figure(figsize=(11, 7), dpi=100)
    # plotting the pressure field as a contour
    pyplot.contourf(Y, X, p, alpha=0.5, cmap=cm.viridis)
    clb = pyplot.colorbar()
    clb.set_label('Pressure')
    # plotting the pressure field outlines
    pyplot.contour(Y, X, p, cmap=cm.viridis)
    # plotting velocity field
    pyplot.quiver(Y[::2, ::2], X[::2, ::2], u[::2, ::2], v[::2, ::2])
    pyplot.xlabel('X')
    pyplot.ylabel('Y')

    rect1 = pyplot.Rectangle((29*dx, 29*dy), 10*dx, 10*dy, color = 'purple')
    rect2 = pyplot.Rectangle((9*dx, 19*dy), 10*dx, 10*dy, color = 'purple')
    rect3 = pyplot.Rectangle((9 * dx, 9 * dy), 5 * dx, 5 * dy, color='purple')

    pyplot.gca().add_patch(rect1)
    pyplot.gca().add_patch(rect2)
    pyplot.gca().add_patch(rect3)

if(plt == 4):
    a = numpy.zeros((nx, ny))
    a[:, :] = (u[:, :]*u[:, :] + v[:, :]*v[:, :])**(1/2)*.5
    pyplot.imshow(a)

if (plt == 5):
    a = numpy.zeros((nx, ny))
    a[:, :] = p[:, :]
    pyplot.imshow(a)

pyplot.show()