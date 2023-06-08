import numpy as np
from matplotlib import pyplot

Nx = 240
L = 6
x_vec = np.linspace(0, L, Nx)
dx = L/(Nx-1)

Nt = 10000
dt = 0.001

h = np.zeros(Nx) + 1
h = np.exp( -((x_vec-(L/2))**2)/.4 ) + 1
hu = np.zeros(Nx)

hb = np.exp( -((x_vec-(L/2))**2)/.4 ) + 1
hub = np.zeros(Nx)

g = 9.8

def difx(s, ds):
    return (  (s[3:] - s[:-3])/(4*ds)  )
def difxx(s, ds):
    return (  (s[2:] - 2*s[1:-1] + s[:-2])/(ds**2)  )

for t in range(0, Nt):
    h0 = h
    hu0 = hu

    h0[0] = h0[-1]
    h0[-1] = h0[-2]
    hu0[0] = -hu0[-1]
    hu0[-1] = -hu0[-2]

    u0 = hu0/h0
    omega0 = (u0**2)*h0 + (1/2)*g*(h0**2)

    # Second order Taylor Series
    h[2:-2] = h0[2:-2] - (dt)*(difx(hu0, dx)) + ((dt**2)/2)*(difxx(omega0, dx))
    hu[2:-2] = hu0[2:-2] - (dt)*(difx(omega0, dx)) + ((dt**2)/2)*(
        -difx(
        -2*u0[2:-2]*difx(omega0, dx) + (u0[2:-2]**2)*difx(hu0, dx) - g*h[2:-2]*difx(hu0, dx)
        , dx)
    )

    # First order Taylor Series
    hb[1:-1] = h0[1:-1] - dt * (difx(hu0, dx))
    hub[1:-1] = hu0[1:-1] - dt * (difx(omega0, dx))

    if(t%1 == 0):
        pyplot.ylim([1, 3])
        pyplot.plot(x_vec, h)
        pyplot.plot(x_vec, hb)
        pyplot.pause(.001)
        pyplot.cla()
pyplot.show()


 