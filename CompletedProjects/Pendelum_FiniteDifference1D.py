import numpy
from matplotlib import pyplot

nt = 200
t_max = 200
t_vec = numpy.linspace(0, t_max, nt)
dt = t_vec[2] - t_vec[1]

nu = .02
g = 9.8
l = 2*2

u = numpy.zeros(nt)

u[0] = 4
u[1] = u[0] + 1

for t in range(1, len(t_vec)-1):
    u[t+1] = dt**2 * ((g/l * numpy.sin(numpy.deg2rad(u[t]))) - (nu * ((u[t] - u[t-1])/dt))) + 2 * u[t] - u[t-1]

    x1 = 0
    x2 = l*numpy.sin(numpy.deg2rad(u[t]))
    y1 = 0
    y2 = l*numpy.cos(numpy.deg2rad(u[t]))

    '''pyplot.plot([x1, x2], [y1, y2])
    pyplot.plot(x2, y2, 'bo')
    pyplot.xlim([-l-1, l+1])
    pyplot.ylim([-l-1, l+1])
    pyplot.pause(.01)
    pyplot.close()'''

pyplot.plot(t_vec, u)
pyplot.show()

