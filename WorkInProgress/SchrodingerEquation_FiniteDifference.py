import numpy as np
from matplotlib import pyplot

Nx = 301
L = 1
Nt = 100000
dx = L/(Nx-1)
dt = 1e-7

x_vec = np.linspace(0, L, Nx)

psi_0 = np.exp(-((x_vec-(L/4))/.3)**2)*2
psi = np.zeros(len(x_vec), dtype=complex)
V = -1e4*np.exp(-(x_vec-(L/2))**2 / (2*(L/20)**2))
V = 1e4*((x_vec-(L/2))/35)**2 - 1.5

for t in range(0, Nt):
    for x in range(1, len(x_vec)-1):
        omega = psi_0[x+1] - 2*psi_0[x] + psi_0[x-1]
        alpha = dt / (dx**2)
        psi[x] = psi_0[x] + 1j/2*alpha*omega - 1j*dt*V[x]*psi_0[x]

    normal = np.sum(np.absolute(psi)**2)*dx
    psi = psi/normal
    psi_0 = psi

    if(t%100==0):
        pyplot.ylim([-2, 6])
        pyplot.xlabel("x/L")
        pyplot.ylabel("|Ψ|^2")
        pyplot.plot(x_vec, np.absolute(psi)**2)
        pyplot.plot(x_vec, V)
        pyplot.pause(.001)
        pyplot.cla()

pyplot.show()








