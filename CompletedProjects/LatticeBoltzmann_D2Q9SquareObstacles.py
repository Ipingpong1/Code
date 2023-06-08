import matplotlib.pyplot as plt
import numpy as np

# Lattice Boltzmann simulation of water passing through a series of boxes in a closed space, includes velocity tracers
# 2022

plot_every = 10
plot_start = 10

def distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

def main():
    Nx = 400 # resolution x-dir
    Ny = 100  # resolution y-dir
    tau = 0.53  # collision timescale
    Nt = 30000  # number of timesteps

    tracer_amount = 25
    tracer_x_initial = 4
    tracer_y_initial = 85

    tracerx = np.zeros(tracer_amount)
    tracery = np.zeros(tracer_amount)
    for i in range(0, tracer_amount):
        tracerx[i] = tracer_x_initial
        tracery[i] = tracer_y_initial-i

    # Lattice speeds / weights
    NL = 9
    cxs = np.array([0, 0, 1, 1, 1, 0, -1, -1, -1])
    cys = np.array([0, 1, 1, 0, -1, -1, -1, 0, 1])
    weights = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])  # sums to 1

    # Initial Conditions
    F = np.ones((Ny, Nx, NL)) # + .01 * np.random.randn(Ny, Nx, NL)
    F[13:86, 0:8, 3] = 2

    # Cylinder boundary
    cylinder = np.full((Ny, Nx), False)

    for y in range(0, Ny):
        for x in range(0, Nx):
            # if (distance(Nx // 4, Ny // 2, x, y) < 13):
            #    cylinder[y][x] = True
            temp = False

            if (y > 0 and y < 10):
                cylinder[y][x] = True
                temp = True

            if (y > 90 and y < 100):
                cylinder[y][x] = True
                temp = True

            if(x>40 and x <100 and y>Ny//1.8):
                cylinder[y][x] = True
                temp = True

            if (x > 140 and x < 200 and y < Ny // 1.8):
                cylinder[y][x] = True
                temp = True

            if (x > 240 and x < 300 and y > Ny // 1.8):
                cylinder[y][x] = True
                temp = True

            if temp == False and y%4 == 0 and x%4 == 0:
                # tracerx = np.append(tracerx, [x])
                # tracery = np.append(tracery, [y])
                pass

    # Simulation Main Loop
    for it in range(Nt):
        print(it)

        F[:, -1, [6, 7, 8]] = F[:, -2, [6, 7, 8]]  # right boundary
        F[:, 0, [2, 3, 4]] = F[:, 1, [2, 3, 4]]  # left boundary

        F[-1, :, [4, 5, 6]] = F[-2, :, [4, 5, 6]]  # up boundary
        F[0, :, [8, 1, 2]] = F[1, :, [8, 1, 2]]  # down boundary

        # Drift
        for i, cx, cy in zip(range(NL), cxs, cys):
            F[:, :, i] = np.roll(F[:, :, i], cx, axis=1)
            F[:, :, i] = np.roll(F[:, :, i], cy, axis=0)

        # Calculate fluid variables
        rho = np.sum(F, 2)
        ux = np.sum(F * cxs, 2) / rho
        uy = np.sum(F * cys, 2) / rho

        # Set reflective boundaries
        bndryF = F[cylinder, :]
        bndryF = bndryF[:, [0, 5, 6, 7, 8, 1, 2, 3, 4]]
        # Apply boundary
        F[cylinder, :] = bndryF
        ux[cylinder] = 0
        uy[cylinder] = 0

        # Apply Collision
        Feq = np.zeros(F.shape)
        for i, cx, cy, w in zip(range(NL), cxs, cys, weights):
            Feq[:, :, i] = rho * w * (
                    1 + 3 * (cx * ux + cy * uy) + 9 * (cx * ux + cy * uy) ** 2 / 2 - 3 * (ux ** 2 + uy ** 2) / 2)

        F = F + -(1.0 / tau) * (F - Feq)

        # ------ TRACER STUFF --------- #
        if it > plot_start:
            for i in range(0, len(tracerx)):
                if (tracerx[i] >= Nx-1):
                    tracerx[i] = 0
                    tracery[i] = 0

                tracerx[i] += ux[int(tracery[i]), int(tracerx[i])]
                tracery[i] += uy[int(tracery[i]), int(tracerx[i])]
            if (it%75 == 0):
                for i in range(0, tracer_amount):
                    tracerx = np.append(tracerx, [tracer_x_initial], axis=0)
                    tracery = np.append(tracery, [(tracer_y_initial - i*3)], axis=0)
        # ----------------------------- #

        if it % plot_every == 0 and it >= plot_start:

            if it == plot_start:
                plt.style.use('dark_background')
                fig, axs = plt.subplots(1)

            fig.suptitle('Lattice-Bolztmann Flow')

            axs.imshow(np.sqrt(uy ** 2 + ux ** 2), cmap=plt.get_cmap('magma'))
            axs.set_title("Velocity Magnitude")

            trace = axs.plot(tracerx, tracery, "wo", markersize="1")

            plt.pause(.001)

            axs.clear()
            plt.cla()

if __name__ == "__main__":
    main()
