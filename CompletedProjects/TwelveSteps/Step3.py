import numpy as nm
from matplotlib import pyplot

v = .3

total_x_domain = 2
x_steps = 40
x_vec = nm.linspace(0, total_x_domain, x_steps)
dx = x_vec[2]-x_vec[1]

total_time = 40
sigma = .2
dt = sigma * dx**2 / v
t_vec = nm.linspace(0, total_time, int(total_time/dt))

u_mat = nm.ones([len(t_vec), len(x_vec)])

for t in range(0, len(t_vec)):
    for x in range(0, len(x_vec)):
        if(x>(len(x_vec)/4) and x<(len(x_vec)/2)):
            # u_mat[t][x] = 2
            pass

        u_mat[t][0] = 200
        u_mat[t][-1] = 200

for t in range(0, len(t_vec)-1):
    print(str(t/len(t_vec) * 100) + '%')
    for x in range(1, len(x_vec)-1):
        u_mat[t+1][x] = v*dt/(dx**2) * (u_mat[t][x-1] - 2*u_mat[t][x] + u_mat[t][x+1]) + u_mat[t][x]

print('100%')

for t in range(0, len(t_vec), 10):
    unew = nm.zeros([1, len(u_mat[t])])
    unew[0][:] = u_mat[t]

    pyplot.figure(figsize=(10, 3))
    pyplot.imshow(unew, vmin=0, vmax=200)
    colorbar = pyplot.colorbar()
    colorbar.set_label("Temperature (C˚)")
    pyplot.pause(.0001)
    colorbar.remove()
    pyplot.close()
