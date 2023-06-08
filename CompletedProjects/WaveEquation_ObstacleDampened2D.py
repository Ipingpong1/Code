import numpy
from matplotlib import pyplot
from matplotlib import cm

# Finite difference solution for the wave equation in one dimension
# 2021

x_max = 2
nx = 300
x_vec = numpy.linspace(0, x_max, nx)
dx = x_vec[2] - x_vec[1]

y_max = 2
ny = 300
y_vec = numpy.linspace(0, y_max, ny)
dy = y_vec[2] - y_vec[1]

cfl = .2
c = 1
nu = .2

dt = cfl*dx/c
t_max = 50
t_vec = numpy.linspace(0, t_max, int(t_max/dt))

u = numpy.zeros([len(t_vec), len(x_vec), len(y_vec)])

for t in range(1, len(t_vec)-1):
    print(t/len(t_vec))

    u[t, 20:40, 20:40] = 0

    if(t<10):
        u[t][50][50] = numpy.sin(1/t)

    u[t+1, 1:-1, 1:-1] = dt**2 * ( (u[t, 2:, 1:-1] - 2*u[t, 1:-1, 1:-1] + u[t, 0:-2, 1:-1])/(dx**2) + (u[t, 1:-1, 2:] - 2*u[t, 1:-1, 1:-1] + u[t, 1:-1, 0:-2])/(dy**2) -
                                  nu * ((u[t, 1:-1, 1:-1] - u[t-1, 1:-1, 1:-1])/dt) ) + \
                         2*u[t, 1:-1, 1:-1] - u[t-1, 1:-1, 1:-1]


min = u.min()
print(min)
max = u.max()
print(max)
for t in range(0, len(t_vec)):
    pyplot.imshow(u[t], vmin=-.05, vmax=.05)
    pyplot.pause(.0001)
    pyplot.close()

pyplot.show()
