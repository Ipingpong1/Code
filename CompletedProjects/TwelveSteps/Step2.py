import numpy as nm
from matplotlib import pyplot

c = 1

total_time = 30
time_steps = 200

t_vec = nm.linspace(0, total_time, time_steps)
dt = t_vec[2] - t_vec[1]

total_domain_x = 40
x_steps = 100

x_vec = nm.linspace(0, total_domain_x, x_steps)
dx = x_vec[2] - x_vec[1]

u_mat = nm.ones([len(t_vec), len(x_vec)])

for t in range(0, len(t_vec)):
    for x in range(0, len(x_vec)):
        if x > (len(x_vec) / 4) and (x < (len(x_vec) / 2)):
            u_mat[t][x] = 2

for t in range(0, len(t_vec)-1):
    for x in range(1, len(x_vec)-1):
        u_mat[t+1][x] = u_mat[t][x] - ((u_mat[t][x]/1.15)*dt/dx) * (u_mat[t][x]-u_mat[t][x-1])

for t in range(0, len(t_vec)):
    pyplot.plot(x_vec, u_mat[t])
    pyplot.pause(.02)
    pyplot.clf()

