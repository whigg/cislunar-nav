# library imports
import os

# local imports
from libs.Batch import *
from libs.Dynamics import *


if __name__ == "__main__":
    # Data and dimensions
    sats = parseGmatData("data/LSIC/8sat_Literature.txt", gmatReport=True)
    n = sats[0].end
    l = 3; m = 3  
        
    # Constants
    rad = 1737400                           # m, radius of moon
    w = 2*np.pi / (27.3217 * 24*60*60)      # rad/s, rotation rate of moon

    fig = plt.figure()
    ax = plt.axes()
    iter = 1

    for i in range(iter):
        x0 = np.array([rad*np.sin(np.pi/12), 0, -rad*np.cos(np.pi/12)])
        xstar = np.array([0, 0, 0])
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
            visible, _ = findVisibleSats(x_true[:,j] / 1000, sats, j, elev=0)
            H = np.diag(computeDOP(x_true[:,j] / 1000, visible, j))
            DOP[:,:,j] = np.diag(H[0:3])

            # Compute measurement
            r = R(DOP[:,:,j], t[j])
            # catch inf or negative vals
            r[(r == float('inf')) | (r < 0)] = sys.float_info.max
            r[r == 0] = 0
            y[:,j] = Y(x_true[:,j], r)


        #np.savetxt('test/true.csv', x_true, delimiter=',')

        # run batch filter
        with Batch(t, xstar, x_true, lambda t: statPhi(statA(w), t), y, G, Ht_3x3,
                lambda z: R(DOP[:,:,np.where(t == z)[0][0]], z)) as batch:
            
            batch.evaluate()
            mc, stat = batch.plot(ax, std=True, semilog=False)
            # update initial guess
            print(batch.x[:,-1])
            #np.savetxt('test/est.csv', batch.x, delimiter=',')

            # print best 3-sigma error
            var = np.diag(batch.P[:,:,j])
            print(3 * np.sqrt(var[0] + var[1] + var[2]))


    ax.grid()
    ax.set_ylim(bottom=0, top=25)
    ax.set_xlim(left=0, right=24)   # bound to actual limits
    ax.set_xlabel("Time (hrs)")
    ax.set_ylabel("Error (m)")
    ax.set_title(f"3σ RMS Position Uncertainty")
    ax.legend(handles=[mc, stat])
    plt.show()