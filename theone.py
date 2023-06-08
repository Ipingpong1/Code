import numpy as np
from matplotlib import pyplot

Lx = 4
Ly = Lx
Nx = 80
Ny = Nx

x_vec = np.linspace(0, Lx, Nx)
y_vec = np.linspace(0, Ly, Ny)

dx = x_vec[1] - x_vec[0]
dy = y_vec[1] - y_vec[0]

Nt = 10000
dt = .001

g = 9.8

h = np.zeros([Nx, Ny])
for x in range(0, len(x_vec)):
    for y in range(0, len(y_vec)):
        h[x, y] = np.exp(- ((x_vec[x]-(Lx/2))**2 + (y_vec[y]-(Ly/2))**2 ) )

hu = np.zeros([Nx, Ny])
hv = np.zeros([Nx, Ny])

h_star_p = np.zeros([Nx, Ny])
hu_star_p = np.zeros([Nx, Ny])
hv_star_p = np.zeros([Nx, Ny])
h_bar_p = np.zeros([Nx, Ny])
hu_bar_p = np.zeros([Nx, Ny])
hv_bar_p = np.zeros([Nx, Ny])

h_star_m = np.zeros([Nx, Ny])
hu_star_m = np.zeros([Nx, Ny])
hv_star_m = np.zeros([Nx, Ny])
h_bar_m = np.zeros([Nx, Ny])
hu_bar_m = np.zeros([Nx, Ny])
hv_bar_m = np.zeros([Nx, Ny])

for t in range(0, Nt):
    h0 = h
    hu0 = hu
    hv0  = hv
    u0 = hu0/h0
    v0 = hv0/h0

    omega_x = h0*(u0**2) + (1/2)*g*(h0**2)
    omega_y = h0*(v0**2) + (1/2)*g*(h0**2)

    h_star_p[1:-1, 1:-1] = (1/2)*(h0[2:, 1:-1] + h0[1:-1, 1:-1]) - dt/(2*dx)*(hu0[2:, 1:-1] - hu0[1:-1, 1:-1])
    hu_star_p[1:-1, 1:-1] = (1/2)*(hu0[2:, 1:-1] + hu0[1:-1, 1:-1]) - dt/(2*dx)*(omega_x[2:, 1:-1] - omega_x[1:-1, 1:-1])
    hv_star_p[1:-1, 1:-1] = (1/2)*(hv0[2:, 1:-1] + hv0[1:-1, 1:-1]) - dt/(2*dx)*(hu0[2:, 1:-1]*v0[2:, 1:-1] - hu0[1:-1, 1:-1]*v0[1:-1, 1:-1])

    h_star_m[1:-1, 1:-1] = (1/2)*(h0[1:-1, 1:-1] + h0[:-2, 1:-1]) - dt/(2*dx)*(hu0[1:-1, 1:-1] - hu0[:-2, 1:-1])
    hu_star_m[1:-1, 1:-1] = (1/2)*(hu0[1:-1, 1:-1] + hu0[:-2, 1:-1]) - dt/(2*dx)*(omega_x[1:-1, 1:-1] - omega_x[:-2, 1:-1])
    hv_star_m[1:-1, 1:-1] = (1/2)*(hv0[1:-1, 1:-1] + hv0[:-2, 1:-1]) - dt/(2*dx)*(hu0[1:-1, 1:-1]*v0[1:-1, 1:-1] - hu0[:-2, 1:-1]*v0[:-2, 1:-1])

    h_bar_p[1:-1, 1:-1] = (1/2)*(h0[1:-1, 2:] + h0[1:-1, 1:-1]) - dt/(2*dx)*(hv0[1:-1, 2:] - hv0[1:-1, 1:-1])
    hu_bar_p[1:-1, 1:-1] = (1/2)*(hu0[1:-1, 2:] + hu0[1:-1, 1:-1]) - dt/(2*dx)*(omega_y[1:-1, 2:] - omega_y[1:-1, 1:-1])
    hv_bar_p[1:-1, 1:-1] = (1/2)*(hv0[1:-1, 2:] + hv0[1:-1, 1:-1]) - dt/(2*dx)*(hu0[1:-1, 2:]*v0[1:-1, 2:] - hu0[1:-1, 1:-1]*v0[1:-1, 1:-1])
    
    h_bar_m[1:-1, 1:-1] = (1/2)*(h0[1:-1, 1:-1] + h0[1:-1, :-2]) - dt/(2*dx)*(hv0[1:-1, 1:-1] - hv0[1:-1, :-2])
    hu_bar_m[1:-1, 1:-1] = (1/2)*(hu0[1:-1, 1:-1] + hu0[1:-1, :-2]) - dt/(2*dx)*(omega_y[1:-1, 1:-1] - omega_y[1:-1, :-2])
    hv_bar_m[1:-1, 1:-1] = (1/2)*(hv0[1:-1, 1:-1] + hv0[1:-1, :-2]) - dt/(2*dx)*(hu0[1:-1, 1:-1]*v0[1:-1, 1:-1] - hu0[1:-1, :-2]*v0[1:-1, :-2] )

    o_star_p_x = ((hu_star_p/h_star_p)**2)*h_star_p + (1/2)*g*(h_star_p**2)
    o_star_m_x = ((hu_star_m/h_star_m)**2)*h_star_m + (1/2)*g*(h_star_m**2)
    o_bar_p_y = ((hv_bar_p/h_bar_p)**2)*h_bar_p + (1/2)*g*(h_bar_p**2)
    o_bar_m_y = ((hv_bar_m/h_bar_m)**2)*h_bar_m + (1/2)*g*(h_bar_m**2)

    h[1:-1, 1:-1] = h0[1:-1, 1:-1] - (dt/dx)*(hu_star_p[1:-1, 1:-1] - hu_star_m[1:-1, 1:-1]) - (dt/dy)*(hv_bar_p[1:-1, 1:-1] - hv_bar_m[1:-1, 1:-1])

    hu[1:-1, 1:-1] = hu0[1:-1, 1:-1] - (dt/dx)*(o_star_p_x[1:-1, 1:-1] - o_star_m_x[1:-1, 1:-1]) - (dt/dy)*(hu_bar_p[1:-1, 1:-1]*(hv_bar_p[1:-1, 1:-1]/h_bar_p[1:-1, 1:-1]) - hu_bar_m[1:-1, 1:-1]*(hv_bar_m[1:-1, 1:-1]/h_bar_m[1:-1, 1:-1]))
    
    hv[1:-1, 1:-1] = hv0[1:-1, 1:-1] - (dt/dx)*(hu_star_p[1:-1, 1:-1]*(hv_star_p[1:-1, 1:-1]/h_star_p[1:-1, 1:-1]) - hu_star_m[1:-1, 1:-1]*(hv_star_m[1:-1, 1:-1]/h_star_m[1:-1, 1:-1]) ) 
    - (dt/dy)*(o_bar_p_y[1:-1, 1:-1] - o_bar_m_y[1:-1, 1:-1])

    if(t%100==0):
        pyplot.imshow(h, vmax=1, vmin=0)
        pyplot.pause(.1)
        pyplot.cla()

pyplot.show()
    






