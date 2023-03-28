# local imports
from filters.LinearKalman import *
from filters.Dynamics import *


if __name__ == "__main__":
    # Data and dimensions
    sats = parseGmatData("data/const_eph/7_MoonOrb.txt", gmatReport=True)
    n = sats[0].end
    l = 3; m = 3; n = 126    
        
    # Constants
    rad = 1737400                           # m, radius of moon
    w = 2*np.pi / (27.3217 * 24*60*60)      # rad/s, rotation rate of moon

    fig = plt.figure()
    ax = plt.axes()
    iter = 10

    for i in range(iter):
        x0 = np.array([np.sin(np.pi/12)*rad, 0., -np.cos(np.pi/12)*rad])
        vx0 = np.array([15**2, 15**2, 15**2])
        xstar = np.array([0, 0, -rad])
        x_true = np.zeros((m,n))
        x_nom  = np.zeros((m,n))
        y = np.zeros((l,n))

        # Dilution of precision
        DOP = np.zeros((l,l,n))
        t = np.linspace(0, 24*60*60, n)    # time steps (in seconds)

        for j in range(n):
            # Compute true trajectory
            x_true[:,j] = statPhi_true(w, t[j]) @ x0

            # Compute DOP at each step
            visible, _ = findVisibleSats(x_true[:,j], sats, j, elev=10)
            H = np.diag(computeDOP(np.array(x0 / 1000), visible, j))
            DOP[:,:,j] = np.diag(H[0:3])

            # Compute measurement
            r = R(DOP[:,:,j], t[j])
            # catch inf or negative vals
            r[(r == float('inf')) | (r < 0)] = sys.float_info.max
            r[r == 0] = 0
            y[:,j] = Y(x_true[:,j], r)

        #np.savetxt('test/true.csv', x_true, delimiter=',')

        # run linear kalman filter
        with LinearKalman(t, x0, np.diag(vx0), x_true, statA, statB,
                        linU, y, Ht_3x3(0), lambda z: linQ(statB(z), (3e-1)**2),
                        lambda z: R(DOP[:,:,np.where(t == z)[0][0]], z)) as dyn:
            
            dyn.evaluate()
            mc, stat = dyn.plot(ax, last=True if i == iter - 1 else False, batch=False, semilog=True)
            # update initial guess
            print(dyn.x[:,-1])
            #np.savetxt('test/est.csv', dyn.x, delimiter=',')


    ax.grid()
    ax.set_ylim(bottom=1, top=1e9)
    ax.set_xlim(left=0, right=24)   # bound to actual limits
    ax.set_xlabel("Time (hrs)")
    ax.set_ylabel("Error (m)")
    ax.set_title(f"RMS Position Uncertainty (4 satellites)")
    ax.legend(handles=[mc, stat])
    plt.show()