import numpy
from matplotlib import pyplot

c = 3

t_total = 4
nt = 2000
t_vec = numpy.linspace(0, t_total, nt)
dt = t_vec[2] - t_vec[1]

x_total = 3
nx = 80
x_vec = numpy.linspace(0, x_total, nx)
dx = x_vec[2] - x_vec[1]

u = numpy.zeros([len(t_vec), len(x_vec)])
# u[:, int(nx/4):int(nx/2)] = .5*(u[:, int(nx/4):int(nx/2)]-.5)*(u[:, int(nx/4):int(nx/2)]-1)
u[0, int(nx/4):int(nx/2)] = 2

dudx = numpy.zeros([len(t_vec), len(x_vec)])

dudt = numpy.zeros([len(t_vec)])

for t in range(0, nt-1):
    u[t+1, 1:-1] = - (c * dt / dx) * (
        u[t, 1:-1] - u[t, :-2]
    ) + u[t, 1:-1]

    dudx[t, 1:-1] = (u[t, 1:-1] - u[t, :-2])/dx

figure, axis = pyplot.subplots(2)
for t in range(0, nt):

    axis[0].plot(x_vec, dudx[t])
    axis[0].set_ylim([-20, 20])

    axis[1].plot(x_vec, u[t])
    axis[1].set_ylim([-.1, 3])

    pyplot.pause(.001)
    axis[0].clear()
    axis[1].clear()

    if(u[t][int(2/dx)]!=0):
        print(t)

pyplot.show()






