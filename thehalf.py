import numpy as np
from matplotlib import pyplot

L = 2
Nx = 80
x_vec = np.linspace(0, L, Nx)
dx = x_vec[2] - x_vec[1]

Nt = 10000
dt = .005

g = 9.8

h = np.zeros(Nx)
h = np.exp( -((x_vec-(L/4))**2)/.4 ) + 1
uh = np.zeros(Nx)

hm = np.zeros(Nx)
uhm = np.zeros(Nx)

for t in range(0, Nt):
    h0 = h
    uh0 = uh
    u0 = uh0/h0

    hm[0:-1] = .5 * (h[0:-1]+h[1:]) - \
            dt/(2*dx) * (uh[1:]-uh[0:-1])

    uhm[0:-1] = .5 * (uh[0:-1]+uh[1:]) - \
             dt / (2 * dx) * ( ((uh[1:]/h[1:])**2*h[1:]+.5*g*h[1:]**2)  -  ((uh[0:-1]/h[0:-1])**2*h[0:-1]+.5*g*h[0:-1]**2)  )

    #  Take a full time step, evaluating the derivative at the half time step,
    #  to estimate the solution at the NX-2 nodes.
    h[1:-1] = h[1:-1] - \
              (dt / dx) * (uhm[1:-1] - uhm[:-2])

    uh[1:-1] = uh[1:-1] - \
               ( dt / dx ) * ( ( ((uhm[1:-1]/hm[1:-1])**2)*hm[1:-1]+.5*g*(hm[1:-1]**2) ) -
                               ( ((uhm[:-2]/hm[:-2])**2)*hm[:-2]+.5*g*(hm[:-2]**2) ) 
                             )
    
    h[0] = h[-1]
    h[-1] = h[-2]
    uh[0] = -uh[-1]
    uh[-1] = -uh[-2]

    pyplot.ylim([0, 3])
    pyplot.xlim([0, L])
    pyplot.plot(x_vec, h)
    pyplot.pause(.01)
    pyplot.clf()

pyplot.show()


