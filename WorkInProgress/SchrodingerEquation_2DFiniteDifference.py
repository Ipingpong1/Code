import numpy as np
from matplotlib import pyplot

Nt = 100000

Nx = 301
Ny = 301
Lx = 1
Ly = 1

x_vec = np.linspace(0, Lx, Nx)
y_vec = np.linspace(0, Ly, Ny)
dx = 1/(Nx-1)
dy = 1/(Ny-1)
dt = 1e-7

# initial conditions
psi = np.zeros([len(x_vec), len(y_vec)], dtype=complex)
psi_0 = np.zeros([len(x_vec), len(y_vec)], dtype=complex)
for x in range(0, len(x_vec)):
    for y in range(0, len(x_vec)):
        psi_0[x, y] = 1.2*np.exp((-(x_vec[x]-1/2)**2 - (y_vec[y]-1/2)**2)/(.1234**2)) # Gaussian Wavefunction

V = np.zeros([len(x_vec), len(y_vec)])
for x in range(0, len(x_vec)):
    for y in range(0, len(y_vec)):
        V[x, y] = -1e5*np.exp(-(x_vec[x]-Lx/2)**2 - (y_vec[y]-Ly/2)**2) # Negative Gaussian Potential Well

for t in range(1, Nt):
    omega_x = (psi_0[2:, 1:-1] - 2*psi_0[1:-1, 1:-1] + psi_0[:-2, 1:-1])/(dx**2)
    omega_y = (psi_0[1:-1, 2:] - 2*psi_0[1:-1, 1:-1] + psi_0[1:-1, :-2])/(dy**2)
    omega = omega_x + omega_y

    psi[1:-1, 1:-1] = psi_0[1:-1, 1:-1] + 1j*dt*omega*.5 - 1j*dt*V[1:-1, 1:-1]*psi_0[1:-1, 1:-1]
    psi_0 = psi

    if(t%100==0):
        pyplot.imshow(np.absolute(psi)**2, vmax=.5)
        pyplot.pause(.0001)
        pyplot.cla()

pyplot.show()


