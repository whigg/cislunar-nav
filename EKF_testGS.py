# local imports
from filters.ExtendedKalman import *
from filters.Dynamics import *


if __name__ == "__main__":
    # Data and dimensions
    sats = parseGmatData("data/const_eph/4_MoonOrb.txt", gmatReport=True)
    n = sats[0].end
    l = 3; m = 6

    np.random.seed(69)
        
    # Constants
    rad = 1737.4                            # m, radius of moon
    ang = 15
    w = 2*np.pi / (27.3217 * 24*60*60)      # rad/s, rotation rate of moon
    W = np.array([0,0,w])
    g = 1.625e-3

    fig = plt.figure()
    ax = plt.axes()
    iter = 1

    for i in range(iter):
        x0 = np.array([rad*np.sin(ang*np.pi/180), 0, -rad*np.cos(ang*np.pi/180), 
            0., w*rad*np.sin(ang*np.pi/180), 0.])
        vx0 = np.array([10**2, 10**2, 10**2, 0.01**2, 0.01**2, 0.01**2])
        # vx0 = np.zeros(np.shape(vx0))

        # Compute true trajectory
        t = np.linspace(0, 24*60*60, n)    # time steps (in seconds)
        u = lambda _: np.zeros((3,))
        func = lambda t, x, u: surfDyn(t, x, u, g, rad, W)
        x_true = integrate(lambda t, x: func(t, x, u), t, x0)

        xstar = np.random.normal(loc=x0, scale=np.sqrt(vx0))
        xstar[0:3] = xstar[0:3] / np.linalg.norm(xstar[0:3]) * rad  # place on surface
        y = np.zeros((l,n))

        # Dilution of precision
        DOP = np.zeros((l,l,n))

        for j in range(n):

            # Compute DOP at each step
            visible, _ = findVisibleSats(x_true[0:3,j], sats, j, elev=0)
            H = np.copy(np.diag(computeDOP(np.array(x0[0:3]), visible, j)))
            # catch negative vals
            H[H <= 0] = float('inf')
            DOP[:,:,j] = np.diag(H[0:3])

            # Compute measurement
            r = R(DOP[:,:,j], t[j])/1e6     # adjust to km
            r[r == float('inf')] = sys.float_info.max   # catch 'inf'
            y[:,j] = Y(x_true[:,j], r)

        np.savetxt('matlab/true.csv', x_true, delimiter=',')

        # run linear kalman filter
        with ExtendedKalman(t, xstar, np.diag(vx0), x_true, func,
                        linU, y, Ht_3x6(0), lambda z: extQ(z, (1e-6)**2),
                        lambda z: R(DOP[:,:,np.where(t == z)[0][0]], z)/1e6) as dyn:
            
            dyn.evaluate()
            mc, stat = dyn.plot(ax, last=True if i == iter - 1 else False, batch=False, semilog=False)
            # update initial guess
            print(dyn.x[:,-1])
            np.savetxt('matlab/est.csv', dyn.x, delimiter=',')


    ax.grid()
    ax.set_ylim(bottom=0, top=0.2)
    ax.set_xlim(left=0, right=24)   # bound to actual limits
    ax.set_xlabel("Time (hrs)")
    ax.set_ylabel("Error (km)")
    ax.set_title(f"RMS Position Uncertainty (4 satellites)")
    ax.legend(handles=[mc, stat])
    plt.show()