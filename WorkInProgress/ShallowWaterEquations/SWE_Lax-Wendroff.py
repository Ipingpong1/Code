from matplotlib import pyplot as plt
import numpy as np

# One Dimension shallow water equation simulation using Lax-Wendroff Scheme
# 2022

nx = 200
nt = 8000
x_length = 1.0
t_length = 1.0
g = 9.8

h = np.zeros(nx) # stores height values
uh = np.zeros(nx) # stores mass velocity values
hm = np.zeros(nx) # stores height values as t + 1/2
uhm = np.zeros(nx) # stores mass velocity values as t + 1/2

h_array = np.zeros([nt, nx]) # stores all height values at all times
uh_array = np.zeros([nt, nx]) # stores all mass velocity values at all times

x_vec = np.linspace(0, x_length, nx) # discretized space domain
t_vec = np.linspace(0, t_length, nt) # discretized time domain

dx = x_vec[2]-x_vec[1] # spacing in space
dt = t_vec[2]-t_vec[1] # spacing in time

# initial conditions
h += 1
h[int(nx / 2 - 5):int(nx / 2 + 5)] = 5

# store the first time step
h_array[0, :] = h[:]
uh_array[0, :] = uh[:]

# main loop
for t in range(0, nt):
    uh*=.9985

    #  Take a half time step, estimating H and UH at the NX-1 spatial midpoints.
    hm[1:-1] = .5 * (h[1:-1]+h[2:]) - \
            dt/(2*dx) * (uh[2:]-uh[1:-1])

    uhm[1:-1] = .5 * (uh[1:-1]+uh[2:]) - \
             dt / (2 * dx) * ( ((uh[2:]/h[2:])**2*h[2:]+.5*g*h[2:]**2)  -  ((uh[1:-1]/h[1:-1])**2*h[1:-1]+.5*g*h[1:-1]**2)  )

    #  Take a full time step, evaluating the derivative at the half time step,
    #  to estimate the solution at the NX-2 nodes.
    h[2:-1] = h[2:-1] - \
              (dt / dx) * (uhm[2:-1] - uhm[1:-2])

    uh[2:-1] = uh[2:-1] - \
               (dt / dx) * ( ((uhm[2:-1]/hm[2:-1])**2*hm[2:-1]+.5*g*hm[2:-1]**2) - ((uhm[1:-2]/hm[1:-2])**2*hm[1:-2]+.5*g*hm[1:-2]**2) )

    # boundary conditions
    uh[:20] = 0
    uh[-20:] = 0

    # add to main array
    h_array[t, :] = h[:]
    uh_array[t, :] = uh[:]

for t in range(0, nt, 80):
    plt.plot(h_array[t])
    plt.xlim([20, nx-20])
    plt.ylim([0, 3.5])
    plt.pause(.001)
    plt.cla()