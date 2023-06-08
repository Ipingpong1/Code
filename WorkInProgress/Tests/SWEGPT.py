import numpy as np
from matplotlib import pyplot as plt

g = 9.8

Lx = 4
Nx = 120
dx = Lx/(Nx-1)

Nt = 10000
dt = .001

x_vec = np.linspace(0, Lx, Nx)
t_vec = np.linspace(0, Nt*dt, Nt)

hu = np.zeros(len(x_vec))
h = np.zeros(len(x_vec)) + 1 + np.exp(-(x_vec-(1))**2/(.4**2)) + np.exp(-(x_vec-(3))**2/(.4**2))
h = np.exp( -((x_vec-(Lx/2))**2)/.4 ) + 1

hu_0 = np.zeros((len(x_vec)))
h_0 = np.zeros((len(x_vec)))
    
def u(hu, h):
    return hu/h
def phi(hu, h, g):
    return h*(u(hu, h)**2) + 1/2*g*(h**2)

nu = 0.000
for t in range(0, len(t_vec)):
    hu_0 = hu
    h_0 = h


    h_0[0] = h_0[-1]
    h_0[-1] = h_0[-2]
    hu_0[0] = -hu_0[-1]
    hu_0[-1] = -hu_0[-2]

    '''
    hu_0[0] = 0
    hu_0[-1] = 0
    h_0[0] = h_0[-1]
    h_0[-1] = h_0[-2]
    u_0 = u(hu_0, h_0)
    '''

    h[1:-1] = -(dt/(2*dx))*(hu_0[2:]-hu_0[:-2]) + h_0[1:-1]
    hu[1:-1] = -(dt/(2*dx))* ( phi(hu_0, h_0, g)[2:]-phi(hu_0, h_0, g)[:-2] - 
    nu*(hu_0[2:] - 2*hu_0[1:-1] + hu_0[:-2])/(dx**2) ) + hu_0[1:-1]

    plt.ylim([0, 2])
    plt.fill_between(x_vec,h)
    plt.plot(x_vec, h)

    plt.pause(.001)
    plt.close()

plt.show()






