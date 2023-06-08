import numpy as np
from matplotlib import pyplot as plt

fig = plt.figure()
ax = plt.axes(projection='3d')

g = 9.8

Lx = 4
Ly = 4

Nx = 40
Ny = 40

x_vec = np.linspace(0, Lx, Nx)
y_vec = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x_vec, y_vec)

dx = x_vec[2] - x_vec[1]
dy = y_vec[2] - y_vec[1]

dt = .00001
Nt = 100000

h = np.zeros([len(x_vec), len(y_vec)]) + 1
for x in range(0, len(x_vec)):
    for y in range(0, len(y_vec)):
        h[x, y] = 2.*np.exp(-((x_vec[x]-(Lx/2))**2+(y_vec[y]-(Lx/2))**2)/(.3))
hu = np.zeros([len(x_vec), len(y_vec)])
hv = np.zeros([len(x_vec), len(y_vec)])

nu = 0.0
for t in range(0, Nt):
    print(t)
    h_0 = h
    hu_0 = hu
    hv_0 = hv

    '''
    hu_0[0, :] = 0
    hu_0[-1, :] = 0
    hu_0[:, 0] = 0
    hu_0[:, -1] = 0

    hv_0[0, :] = 0
    hv_0[-1, :] = 0
    hv_0[:, 0] = 0
    hv_0[:, -1] = 0

    h_0[0, :] = h_0[-1, :]
    h_0[-1, :] = h_0[-2, :]
    h_0[:, 0] = h_0[:, -1]
    h_0[:, -1] = h_0[:, -2]'''

    u = hu_0/h_0
    v = hv_0/h_0
    e = h_0*u*v
    o1 = ( h_0*(u**2) + .5*g*(h_0**2) )
    o2 = ( h_0*(v**2) + .5*g*(h_0**2) )

    h[1:-1, 1:-1] = -dt * ( 
        (hu_0[2:, 1:-1]-hu_0[:-2, 1:-1])/(2*dx) + \
        (hv_0[1:-1, 2:]-hv_0[1:-1, :-2])/(2*dy) 
        ) + h_0[1:-1, 1:-1] 

    hu[1:-1, 1:-1] = -dt * (
        (o1[2:, 1:-1] - o1[:-2, 1:-1])/(2*dx) + \
        (e[1:-1, 2:] - e[1:-1, :-2])/(2*dy) - \
        nu*( ((hu_0[1:-1, 2:] - 2*hu_0[1:-1, 1:-1] + hu_0[1:-1, :-2]) +
              (hu_0[2:, 1:-1] - 2*hu_0[1:-1, 1:-1] + hu_0[:-2, 1:-1])) /(dx**2) ) ) + hu_0[1:-1, 1:-1]

    hv[1:-1, 1:-1] = -dt * (
        (o2[1:-1, 2:] - o2[1:-1, :-2])/(2*dy) + \
        (e[2:, 1:-1] - e[:-2, 1:-1])/(2*dx) - \
        nu*( ((hv_0[1:-1, 2:] - 2*hv_0[1:-1, 1:-1] + hv_0[1:-1, :-2]) +
              (hv_0[2:, 1:-1] - 2*hv_0[1:-1, 1:-1] + hv_0[:-2, 1:-1]))/(dy**2) ) ) + hv_0[1:-1, 1:-1]

    if(t>55000):
        #surf = ax.plot_surface(X, Y, h, color='b', shade=True,
        #                       linewidth=0, antialiased=False)
        ax.plot_surface(X, Y, h, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
        # ax.view_init(elev=45)
        ax.set_zlim(-.0001, 2)
        plt.axis('off')
        plt.pause(.0001)
        plt.cla()

    




