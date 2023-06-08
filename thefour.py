import numpy as np
from matplotlib import pyplot

Nx = 60
Ny = Nx
L = 6

x_vec = np.linspace(0, L, Nx)
y_vec = np.linspace(0, L, Ny)
dx = L/(Nx-1)
dy = L/(Ny-1)

Nt = 10000
dt = 0.0005

h = np.zeros([Nx, Ny])
for x in range(0, len(x_vec)):
    for y in range(0, len(y_vec)):
        h[x, y] =  1 + .4*np.exp(- ((x_vec[x]-(L/2))**2 + (y_vec[y]-(L/2))**2 ) )
hu = np.zeros([Nx, Ny])
hv = np.zeros([Nx, Ny])

g = 9.8

def difx(s):
    return (  (s[1:, :-1] - s[:-1, :-1]) / (dx)  )
def dify(s):
    return (  (s[:-1, 1:] - s[:-1, :-1]) / (dy)  )

for t in range(0, Nt):
    h0 = h
    hu0 = hu
    hv0 = hv

    hu0 *= .999
    hv0 *= .999

    hu0[0, :] = 0 #- hu0[-1, :]
    hu0[-1, :] = 0 #- hu0[-2, :]
    hu0[:, 0] = 0 #- hu0[:, -1]
    hu0[:, -1] = 0 #- hu0[:, -2]

    hv0[0, :] = 0 #- hv0[-1, :]
    hv0[-1, :] = 0 #- hv0[-2, :]
    hv0[:, 0] = 0 #- hv0[:, -1]
    hv0[:, -1] = 0 #- hv0[:, -2]

    h0[0, :] = 1 #h0[-1, :]
    h0[-1, :] = 1 #h0[-2, :]
    h0[:, 0] = 1 #h0[:, -1]
    h0[:, -1] = 1 #h0[:, -2]

    u0 = hu0/h0
    v0 = hv0/h0
    omegax = (u0**2)*h0 + (1/2)*g*(h0**2)
    omegay = (v0**2)*h0 + (1/2)*g*(h0**2)
    huv = hu0*hv0/h0

    h[0:-1, 0:-1]  =  h0[0:-1, 0:-1] - dt * ( difx(hu0) + dify(hv0) )
    hu[0:-1, 0:-1] = hu0[0:-1, 0:-1] - dt * ( difx(omegax) + dify(huv) )
    hv[0:-1, 0:-1] = hv0[0:-1, 0:-1] - dt * ( dify(omegay) + difx(huv) )

    if(t%100 == 0):
        disp = np.zeros([int(Nx/2 - 1), int(Ny/2 - 1)])
        for x in range(0, Nx-1):
            for y in range(0, Ny-1):
                if((x+1)%2 == 0 and (y+1)%2 == 0):
                    disp[int(x/2), int(y/2)] = h[x, y]

        pyplot.imshow(h)
        pyplot.pause(.001)
        pyplot.cla()
pyplot.show()


 