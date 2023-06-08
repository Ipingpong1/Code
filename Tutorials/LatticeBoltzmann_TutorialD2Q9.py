import matplotlib.pyplot as plt
import numpy as np

plot_every = 10
plot_start = 200

def distance(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

def main():
    Nx = 400 # resolution x-dir
    Ny = 100  # resolution y-dir
    tau = 0.53  # collision timescale
    Nt = 30000  # number of timesteps

    # Lattice speeds / weights
    NL = 9
    cxs = np.array([0, 0, 1, 1, 1, 0, -1, -1, -1])
    cys = np.array([0, 1, 1, 0, -1, -1, -1, 0, 1])
    weights = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])  # sums to 1

    # Initial Conditions
    F = np.ones((Ny, Nx, NL)) + .01 * np.random.randn(Ny, Nx, NL)
    F[:, :, 3] = 2.3

    # Cylinder boundary
    cylinder = np.full((Ny, Nx), False)

    for y in range(0, Ny):
        for x in range(0, Nx):
            if (distance(Nx // 4, Ny // 2, x, y) < 13):
                cylinder[y][x] = True

    # Simulation Main Loop
    for it in range(Nt):
        print(it)

        F[:, -1, [6, 7, 8]] = F[:, -2, [6, 7, 8]]  # right boundary
        F[:, 0, [2, 3, 4]] = F[:, 1, [2, 3, 4]]  # left boundary

        # Drift
        for i, cx, cy in zip(range(NL), cxs, cys):
            F[:, :, i] = np.roll(F[:, :, i], cx, axis=1)
            F[:, :, i] = np.roll(F[:, :, i], cy, axis=0)

        # Set reflective boundaries
        bndryF = F[cylinder, :]
        bndryF = bndryF[:, [0, 5, 6, 7, 8, 1, 2, 3, 4]]

        # Calculate fluid variables
        rho = np.sum(F, 2)
        ux = np.sum(F * cxs, 2) / rho
        uy = np.sum(F * cys, 2) / rho

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


        if it % plot_every == 0 and it >= plot_start:

            if it == plot_start:
                plt.style.use('dark_background')
                fig, axs = plt.subplots(1)

            fig.suptitle('Lattice-Bolztmann Flow')

            dfydx = ux[2:, 1:-1] - ux[0:-2, 1:-1]
            dfxdy = uy[1:-1, 2:] - uy[1:-1, 0:-2]
            curlmat = dfydx - dfxdy
            axs.imshow(curlmat, cmap="bwr", vmin=-.1, vmax=.1)
            axs.set_title("Curl Magnitude")

            # axs.imshow(np.sqrt(uy ** 2 + ux ** 2), cmap=plt.get_cmap('magma'))
            # axs.set_title("Velocity Magnitude")

            plt.pause(.001)

            axs.clear()
            plt.cla()

if __name__ == "__main__":
    main()
