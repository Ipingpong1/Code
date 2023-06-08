import numpy
from matplotlib import pyplot
from matplotlib import cm

# Simulation of the 2 dimension wave equation with an square obstacle in its path
# 2021

x_max = 10
x_min = 0
nx = 80
x = numpy.linspace(x_min, x_max, nx)
dx = x[2] - x[1]

y_max = 10
y_min = 0
ny = nx
y = numpy.linspace(y_min, y_max, ny)
dy = y[2] - y[1]

c = 1

CFL = .2

dt = c*CFL*dx
nt = 400

nu = .002

u = numpy.zeros([nt, nx, ny])
u[0, nx // 2, ny // 2] = numpy.sin(0)

for t in range(1, nt-1):
    if(t<100):
       u[t, nx // 2, ny // 2] = numpy.sin(t/10)

    u[t+1, 1:-1, 1:-1] = c * (dt**2/(dx**2)) * ( (u[t, 2:, 1:-1] - 2*u[t, 1:-1, 1:-1] + u[t, :-2, 1:-1]) + (u[t, 1:-1, 2:] - 2*u[t, 1:-1, 1:-1] + u[t, 1:-1, :-2]) -
                                                 nu*(u[t, 1:-1, 1:-1] - u[t-1, 1:-1, 1:-1])/dt ) \
                         + 2*u[t, 1:-1, 1:-1] - u[t-1, 1:-1, 1:-1]

fig = pyplot.figure()
X, Y = numpy.meshgrid(x, y)
ax = fig.add_subplot(111, projection='3d')
for t in range(0, nt):
    surf = ax.plot_surface(X, Y, u[t], color='b', shade=True,
                           linewidth=0, antialiased=False)
    ax.view_init(elev=45)
    ax.set_zlim(-.0001, 2.4)
    pyplot.axis('off')

    pyplot.pause(.0001)
    pyplot.cla()