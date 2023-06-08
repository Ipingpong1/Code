import matplotlib.pyplot as plt
import numpy as np

def shallow_water_1d_test():

    #  Set parameters.
    nx = 81
    nt = 400
    x_length = 1.0
    t_length = 0.2
    g = 9.8

    #  Compute H and UH.
    h_array, uh_array, x, t = shallow_water_1d(nx, nt, x_length, t_length, g)
    h_array = h_array.T

    for t in range(0, nt):
        plt.plot(h_array[t])
        # plt.show()
        plt.xlim([0, nx])
        plt.ylim([0, 3.5])
        plt.pause(.001)
        plt.cla()

    return


def shallow_water_1d(nx, nt, x_length, t_length, g):
    #    This code can be considered a 1D version of Cleve Moler's shallow
    #    water equation solver.
    #
    #    The version of the shallow water equations being solved here is in
    #    conservative form, and omits the Coriolis force.  The state variables
    #    are H (the height) and UH (the mass velocity).
    #
    #    The equations have the form
    #
    #      dH/dt + d UH/dx = 0
    #
    #      d UH/dt + d ( U^2 H + 1/2 g H^2 )/dx = 0
    #
    #    Here U is the ordinary velocity, U = UH/H, and g is the gravitational
    #    acceleration.
    #
    #    The initial conditions are used to specify ( H, UH ) at an equally
    #    spaced set of points, and then the Lax-Wendroff method is used to advance
    #    the solution through a number of equally spaced points in time, with
    #    boundary conditions supplying the first and last spatial values.
    #
    #    Some input values will result in an unstable calculation that
    #    quickly blows up.  This is related to the Courant-Friedrichs-Lewy
    #    condition, which requires that DT be small enough, relative to DX and
    #    the velocity, that information cannot cross an entire cell.
    #
    #    A "reasonable" set of input quantities is
    #
    #      shallow_water_1d ( 41, 100, 1.0, 0.2, 9.8 )
    #
    #  Input:
    #
    #    integer NX, the number of spatial nodes.
    #
    #    integer NT, the number of times steps.
    #
    #    real X_LENGTH, the length of the region.
    #
    #    real T_LENGTH, the time extent.
    #
    #    real G, the gravity constant.  G = 9.8 meters per second^2.
    #
    #  Output:
    #
    #    real H_ARRAY(NX,NT+1), the height for all space and time points.
    #
    #    real UH_ARRAY(NX,NT+1), the mass velocity for all space and time points.
    #
    #    real X(NX), the X coordinates.
    #
    #    real T(NT+1), the T coordinates.

    #  Allocate vectors.
    h = np.zeros(nx)
    uh = np.zeros(nx)
    hm = np.zeros(nx - 1)
    uhm = np.zeros(nx - 1)
    h_array = np.zeros([nx, nt + 1])
    uh_array = np.zeros([nx, nt + 1])

    #  Define the locations of the nodes and time steps and the spacing.
    x = np.linspace(0, x_length, nx)
    t = np.linspace(0, t_length, nt + 1)

    dx = x_length / float(nx - 1)
    dt = t_length / float(nt)

    #  Apply the initial conditions.
    h, uh = initial_conditions(nx, nt, h, uh, x)

    #  Apply the boundary conditions.
    h, uh = boundary_conditions(nx, nt, h, uh, t[0])

    #  Store the first time step into H_ARRAY and UH_ARRAY.
    h_array[0:nx, 0] = h[0:nx]
    uh_array[0:nx, 0] = uh[0:nx]

    #  Take NT more time steps.
    for it in range(1, nt + 1):
        #  Take a half time step, estimating H and UH at the NX-1 spatial midpoints.
        hm[0:nx - 1] = (h[0:nx - 1] + h[1:nx]) / 2.0 \
                       - (dt / 2.0) * (uh[1:nx] - uh[0:nx - 1]) / dx

        uhm[0:nx - 1] = (uh[0:nx - 1] + uh[1:nx]) / 2.0 \
                        - (dt / 2.0) * ( \
                                    uh[1:nx] ** 2 / h[1:nx] + 0.5 * g * h[1:nx] ** 2 \
                                    - uh[0:nx - 1] ** 2 / h[0:nx - 1] - 0.5 * g * h[0:nx - 1] ** 2) / dx

        #  Take a full time step, evaluating the derivative at the half time step,
        #  to estimate the solution at the NX-2 nodes.
        h[1:nx - 1] = h[1:nx - 1] \
                      - dt * (uhm[1:nx - 1] - uhm[0:nx - 2]) / dx

        uh[1:nx - 1] = uh[1:nx - 1] \
                       - dt * ( \
                                   uhm[1:nx - 1] ** 2 / hm[1:nx - 1] + 0.5 * g * hm[1:nx - 1] ** 2 \
                                   - uhm[0:nx - 2] ** 2 / hm[0:nx - 2] - 0.5 * g * hm[0:nx - 2] ** 2) / dx

        #  Update the boundary conditions.
        h, uh = boundary_conditions(nx, nt, h, uh, t[it])

        #  Copy data into the big arrays.
        h_array[0:nx, it] = h[0:nx]
        uh_array[0:nx, it] = uh[0:nx]

    return h_array, uh_array, x, t


def boundary_conditions(nx, nt, h, uh, t):
    #  Input:
    #
    #    integer NX, the number of spatial nodes.
    #
    #    integer NT, the number of times steps.
    #
    #    real H(NX), the height for all space.
    #
    #    real UH(NX), the mass velocity for all space.
    #
    #    real T, the current time.
    #
    #  Output:
    #
    #    real H(NX), the height, with H(1) and H(NX) adjusted for
    #    boundary conditions.
    #
    #    real UH(NX), the mass velocity, with UH(1) and UH(NX)
    #    adjusted for boundary conditions.
    #
    bc = 1

    #  Periodic boundary conditions on H and UH.
    if (bc == 1):
        h[0] = h[nx - 2]
        h[nx - 1] = h[1]
        uh[0] = uh[nx - 2]
        uh[nx - 1] = uh[1]

    #  Free boundary conditions on H and UH.
    elif (bc == 2):
        h[0] = h[1]
        h[nx - 1] = h[nx - 2]
        uh[0] = uh[1]
        uh[nx - 1] = uh[nx - 2]

    #  Reflective boundary conditions on UH, free boundary conditions on H.
    elif (bc == 3):
        h[0] = h[1]
        h[nx - 1] = h[nx - 2]
        uh[0] = - uh[1]
        uh[nx - 1] = - uh[nx - 2]

    return h, uh


def initial_conditions(nx, nt, h, uh, x):

    # h = 2.0 + np.sin(2.0 * np.pi * x)
    h = np.ones(nx)
    h[40] = 2
    uh = np.zeros(nx)

    return h, uh


if (__name__ == '__main__'):
    shallow_water_1d_test()





