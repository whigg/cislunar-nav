# local imports
from filters.LinearKalman import *
from filters.Dynamics import *


if __name__ == "__main__":
    # Data and dimensions
    sats = parseGmatData("gmat/GPS_20161231_24sat_1day.txt", gmatReport=True)
    n = sats[0].end
    l = 3; m = 3; n = 126    
        
    # Constants
    r = 1737400                             # m, radius of moon
    w = 2*np.pi / (27.3217 * 24*60*60)      # rad/s, rotation rate of moon
    # x0 = np.array([np.sin(np.pi/12)*r, 0., -np.cos(np.pi/12)*r, 0, w*np.sin(np.pi/12)*r, 0])
    x0 = np.array([np.sin(np.pi/12)*r, 0., -np.cos(np.pi/12)*r])
    vx0 = np.array([15**2, 15**2, 15**2])
    xstar = np.array([0, 0, -r])
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
        y[:,j] = Y(x_true[:,j], R(DOP[:,:,j], t[j]))

    np.savetxt('test/true.csv', x_true, delimiter=',')

    # run linear kalman filter
    with LinearKalman(t, x0, np.diag(vx0), x_true, statA, statB,
                      linU, y, Ht_3x3(0), lambda z: linQ(statB(z), (3e-1)**2),
                      lambda z: R(DOP[:,:,np.where(t == z)[0][0]], z)
        ) as dyn:
        
        fig = plt.figure()
        ax = plt.axes()
        
        iter = 1
        for i in range(iter):
            dyn.evaluate()
            mc, stat = dyn.plot(ax, last=True if i == iter - 1 else False, batch=False)
            # update initial guess
            print(dyn.x[:,-1])
            #dyn.x0 = dyn.x[:,-1]
        np.savetxt('test/est.csv', dyn.x, delimiter=',')

        ax.grid()
        #ax.set_ylim(bottom=0, top=65)
        ax.set_xlabel("Samples")
        ax.set_ylabel("Error (m)")
        ax.set_title(f"RMS Position Uncertainty (Linear Kalman filter, GPS geometry)")
        ax.legend(handles=[mc, stat])
        plt.show()