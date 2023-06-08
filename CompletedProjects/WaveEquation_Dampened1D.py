import numpy
from matplotlib import pyplot

# Finite difference solution for the wave equation in one dimension
# 2021

x_max = 2
nx = 300
x_vec = numpy.linspace(0, x_max, nx)
dx = x_vec[2] - x_vec[1]

cfl = .2
dt = cfl*dx
t_max = 10
t_vec = numpy.linspace(0, t_max, int(t_max/dt))

nu = 2*numpy.pi/x_max
k = 1

u = numpy.zeros((len(t_vec), len(x_vec)))

for t in range(1, len(t_vec)-1):
    print(t/len(t_vec))
    if(t<20):
        u[t, 151] = 2*numpy.sin(1/t)

    for x in range(1, len(x_vec)-1):
        u[t+1, x] = (k*(dt**2)) * ((u[t, x+1] - 2*u[t, x] + u[t, x-1])/(dx**2) - nu * ((u[t, x] - u[t-1, x])/dt)) - u[t-1, x] + 2*u[t, x]


for t in range(0, len(t_vec)):
    pyplot.xlim([0, x_max])
    pyplot.ylim([-2, 2])
    pyplot.plot(x_vec, u[t])
    pyplot.pause(.001)
    pyplot.cla()

pyplot.show()