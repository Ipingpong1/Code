import numpy as np
from matplotlib import pyplot

fig = pyplot.figure()
ax = pyplot.axes(projection='3d')

Lx = 15
Ly = Lx

Nx = 60
Ny = Nx

x_vec = np.linspace(0, Lx, Nx)
y_vec = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x_vec, y_vec)

dx = x_vec[1] - x_vec[0]
dy = y_vec[1] - y_vec[0]

Nt = 10000
dt = .005

g = 9.8

hu = np.zeros([Nx, Ny])
hv = np.zeros([Nx, Ny])
h = np.zeros([Nx, Ny])
for x in range(0, len(x_vec)):
    for y in range(0, len(y_vec)):
        h[x, y] =  1 + .2*np.exp( (- ((x_vec[x]-(Lx/2))**2 + (y_vec[y]-(Ly/3))**2 ))/.5 )

h_xp = np.copy(h)
hu_xp = np.copy(h)
hv_xp = np.copy(h)

h_xm = np.copy(h)
hu_xm = np.copy(h)
hv_xm = np.copy(h)

h_yp = np.copy(h)
hu_yp = np.copy(h)
hv_yp = np.copy(h)

h_ym = np.copy(h)
hu_ym = np.copy(h)
hv_ym = np.copy(h)

for t in range(Nt):
    h0 = h
    hu0 = hu
    hv0 = hv

    hu0[0, :] = - hu0[-1, :]
    hu0[-1, :] = - hu0[-2, :]
    hu0[:, 0] = - hu0[:, -1]
    hu0[:, -1] = - hu0[:, -2]

    hv0[0, :] = - hv0[-1, :]
    hv0[-1, :] = - hv0[-2, :]
    hv0[:, 0] = - hv0[:, -1]
    hv0[:, -1] = - hv0[:, -2]

    h0[0, :] = h0[-1, :]
    h0[-1, :] = h0[-2, :]
    h0[:, 0] = h0[:, -1]
    h0[:, -1] = h0[:, -2]

    # wx = 6
    # wy = 6
    # hu0[int(Nx/2-wx):int(Nx/2+wx), int(Ny/2-wy):int(Ny/2+wy)] = 0
    # hv0[int(Nx/2-wx):int(Nx/2+wx), int(Ny/2-wy):int(Ny/2+wy)] = 0
    # h0[int(Nx/2-wx):int(Nx/2+wx), int(Ny/2-wy):int(Ny/2+wy)]  = 1

    u0 = hu0/h0
    v0 = hv0/h0

    hu0 *= .998
    hv0 *= .998

    # -------------------- #
    # Calculate half-steps #
    # -------------------- #

    # -------------------- U^{n+1/2}_{i+1/2, j} -------------------- #
    h_xp[0:-1, 0:-1] = (1/2)*(h0[1:, 0:-1] + h0[0:-1, 0:-1]) - \
                    (dt/(2*dx))* \
                    (
                    ( hu0[1:, 0:-1] ) # F^n_{i+1, j} (h)
                    - 
                    ( hu0[0:-1, 0:-1] ) # F^n_{i, j} (h)
                    )
    
    hu_xp[0:-1, 0:-1] = (1/2)*(hu0[1:, 0:-1] + hu0[0:-1, 0:-1]) - \
                    (dt/(2*dx))*\
                    ( 
                    ( (u0[1:, 0:-1]**2) * h0[1:, 0:-1] + (1/2)*g*(h0[1:, 0:-1]**2) ) # F^n_{i+1, j} (hu)
                     -
                    ( (u0[0:-1, 0:-1]**2) * h0[0:-1, 0:-1] + (1/2)*g*(h0[0:-1, 0:-1]**2) ) # F^n_{i, j} (hu)
                    )
    
    hv_xp[0:-1, 0:-1] = (1/2)*(hv0[1:, 0:-1] + hv0[0:-1, 0:-1]) - \
                    (dt/(2*dx))*\
                    (
                    ( hu0[1:, 0:-1]*v0[1:, 0:-1] ) # F^n_{i+1, j} (hv)
                    -
                    ( hu0[0:-1, 0:-1]*v0[0:-1, 0:-1] ) # F^n_{i, j} (hv)
                    )
    

    # -------------------- U^{n+1/2}_{i, j+1/2} -------------------- #
    h_yp[0:-1, 0:-1] = (1/2)*(h0[0:-1, 1:] + h0[0:-1, 0:-1]) - \
                    (dt/(2*dy))*\
                        (
                        hv0[0:-1, 1:] # G^n_{i, j+1} (h)
                        -
                        hv0[0:-1, 0:-1] # G^n_{i, j} (h)
                        )
    
    hu_yp[0:-1, 0:-1] = (1/2)*(hu0[0:-1, 1:] + hu0[0:-1, 0:-1]) - \
                    (dt/(2*dy))*\
                        (
                        ( hu0[0:-1, 1:]*v0[0:-1, 1:] ) # G^n_{i, j+1} (hu)
                        -
                        ( hu0[0:-1, 0:-1]*v0[0:-1, 0:-1] ) # G^n_{i, j} (hu)
                        )
    
    hv_yp[0:-1, 0:-1] = (1/2)*(hv0[0:-1, 1:] + hv0[0:-1, 0:-1]) - \
                    (dt/(2*dy))*\
                        (
                        ( (v0[0:-1, 1:]**2)*h0[0:-1, 1:] + (1/2)*g*(h0[0:-1, 1:]**2) ) # G^n_{i, j+1} (hv)
                        -
                        ( (v0[0:-1, 0:-1]**2)*h0[0:-1, 0:-1] + (1/2)*g*(h0[0:-1, 0:-1]**2) ) # G^n_{i, j} (hv)
                        )
    
    # ------------------------------------------ #
    # Calculate fluxes using half-step variables #
    # ------------------------------------------ #
    #                                            #
    #        /          hu               \ [h]   #
    #    F = | h*(u**2) + (1/2)*g*(h**2) | [hu]  #
    #        \         huv               / [hv]  #
    #                                            #
    #        /          hv               \ [h]   #
    #    G = |          huv              | [hu]  #
    #        \ h*(v**2) + (1/2)*g*(h**2) / [hv]  #
    #                                            #
    # ------------------------------------------ #

    F_h_p = hu_xp
    G_h_p = hv_yp

    F_hu_p = h_xp * ((hu_xp/h_xp)**2) + (1/2)*g*(h_xp**2)
    G_hu_p = hu_yp * (hv_yp/h_yp)

    F_hv_p = hu_xp * (hv_xp/h_xp)
    G_hv_p = h_yp*((hv_yp/h_yp)**2) + (1/2)*g*(h_yp**2)

    # ---------------------------------------- #
    # Calculate final t+dt values using fluxes #
    # ---------------------------------------- #

    h[1:-1, 1:-1] = h0[1:-1, 1:-1] - (dt/dx) * (F_h_p[1:-1, 1:-1] - F_h_p[:-2, 1:-1]) - (dt/dy) * (G_h_p[1:-1, 1:-1] - G_h_p[1:-1, :-2])
    hu[1:-1, 1:-1] = hu0[1:-1, 1:-1] - (dt/dx) * (F_hu_p[1:-1, 1:-1] - F_hu_p[:-2, 1:-1]) - (dt/dy) * (G_hu_p[1:-1, 1:-1] - G_hu_p[1:-1, :-2])
    hv[1:-1, 1:-1] = hv0[1:-1, 1:-1] - (dt/dx) * (F_hv_p[1:-1, 1:-1] - F_hv_p[:-2, 1:-1]) - (dt/dy) * (G_hv_p[1:-1, 1:-1] - G_hv_p[1:-1, :-2])

    # ---------------------------------------------- #
    # Boundary conditions (ref on hu, hv; free on h) #
    # ---------------------------------------------- #

    cw = 5
    loc = 1.5
    ow = 2
    h[:int(Nx/2-ow), int(Nx/loc-cw):int(Nx/loc+cw)] = 1
    h[int(Nx/2+ow):, int(Nx/loc-cw):int(Nx/loc+cw)] = 1

    hu[:int(Nx/2-ow), int(Nx/loc-cw):int(Nx/loc+cw)] = 0
    hu[int(Nx/2+ow):, int(Nx/loc-cw):int(Nx/loc+cw)] = 0

    hv[:int(Nx/2-ow), int(Nx/loc-cw):int(Nx/loc+cw)] = 0
    hv[int(Nx/2+ow):, int(Nx/loc-cw):int(Nx/loc+cw)] = 0

    # ----------------------- #
    # Plot every 20 timesteps #
    # ----------------------- #

    if(t%20==0):
        ax.plot_surface(X, Y, h, cmap="viridis", shade=True, antialiased=False, linewidth=0)
        ax.set_zlim(1, 1.05)
        # pyplot.imshow(h)
        pyplot.pause(.0001)
        pyplot.cla()

pyplot.show()




    
