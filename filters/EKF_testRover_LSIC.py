# library imports
import os
from cycler import cycler

# local imports
from libs.ExtendedKalman import *
from libs.Dynamics import *


if __name__ == "__main__":
    # Data and dimensions
    dir = 'data/LSIC/'              # relative path to data files
    files = []
    for file in os.listdir(dir):    # create list of files in directory
        f = os.path.join(dir, file)

        if os.path.isfile(f) and 'Literature' in f:
            files.append(f)

    np.random.seed(69)
        
    # Constants
    l = 3; m = 6
    rad = 1737.4                            # m, radius of moon
    ang = 15
    w = 2*np.pi / (27.3217 * 24*60*60)      # rad/s, rotation rate of moon
    W = np.array([0,0,w])
    g = 1.625e-3
    std_accel = 1e-10    # root Allan variance

    # plotting
    default_cycler = (cycler(color=['b','g','r','c','m']) + cycler(linestyle=['-','--',':','-.','-']))
    plt.rc('axes', prop_cycle=default_cycler)
    fig = plt.figure()
    ax = plt.axes()
    stdplot = []
    labels = ['4 Sat*', '8 Sat*']
    # labels = ['5 Sat', '6 Sat']

    for i, file in enumerate(files):
        # load data
        sats = parseGmatData(file, gmatReport=True)
        n = sats[0].end

        v0 = 0.000033
        x0 = np.array([rad*np.sin(ang*np.pi/180), 0, -rad*np.cos(ang*np.pi/180), 
            np.cos(ang*np.pi/180)*v0, w*rad*np.sin(ang*np.pi/180), np.sin(ang*np.pi/180)*v0])
        vx0 = np.array([1**2, 1**2, 1**2, 0**2, 0**2, 0**2])
        # vx0 = np.zeros(np.shape(vx0))
      
        # Compute true trajectory
        t = np.linspace(0, 24*60*60, n)    # time steps (in seconds)
        func = lambda t, x, u: surfDyn(t, x, u, g, rad, W)
        randWalk, randWalkErr = getAccelFuncs(t, 10, 1e-7, std_accel)  
        x_true = integrate(lambda t, x: func(t, x, randWalk), t, x0)

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

        # run extended kalman filter
        with ExtendedKalman(t, xstar, np.diag(vx0), x_true, func,
                        randWalkErr, y, Ht_3x6(0), lambda z: linQ(statB_6x3(z), (std_accel * 3)**2),
                        lambda z: R(DOP[:,:,np.where(t == z)[0][0]], z)/1e6) as dyn:
            
            dyn.evaluate()
            dyn.units = 1e3         # unit conversion for plotting
            _, stat = dyn.plot(ax, mc=False, std=True, batch=False, semilog=False)
            stdplot.append(stat)
            # update initial guess
            print(dyn.x[:,-1])
            np.savetxt('matlab/est.csv', dyn.x, delimiter=',')


    ax.grid()
    ax.set_ylim(bottom=0, top=300)
    ax.set_xlim(left=0, right=24)   # bound to actual limits
    ax.set_xlabel("Time (hrs)")
    ax.set_ylabel("Error (m)")
    ax.set_title(f"3σ RMS Position Uncertainty, EKF")
    ax.legend(labels)
    plt.show()