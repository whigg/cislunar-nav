# local imports
from Batch import *
from URE import *

def A(w):
    ''' state dynamics matrix '''
    A11 = np.array([[0, -w, 0],[w, 0, 0],[0, 0, 0]])
    # A1 = np.concatenate((A11, np.zeros((3,3))), axis=1)
    # A2 = np.concatenate((np.zeros((3,3)), A11), axis=1)
    # return np.concatenate((A1,A2), axis=0)
    return A11

def phi(A, t):
    ''' state transition matrix '''
    return expm(A * t)

def phi_true(w, t):
    ''' true state transition matrix '''
    Pd = np.array([[np.cos(w*t), -np.sin(w*t), 0],[np.sin(w*t), np.cos(w*t), 0],[0, 0, 1]])
    # P1 = np.concatenate((Pd,np.zeros((3,3))), axis=1)
    # P2 = np.concatenate((np.zeros((3,3)),Pd), axis=1)
    # return np.concatenate((P1,P2), axis=0)
    return Pd

def Y(X, R):
    ''' true measurement '''
    # X is true state and R is measurement covariance
    return np.random.normal(loc=X[0:3], scale=np.sqrt(np.diag(R)))

def G(X):
    ''' measurement model '''
    return X[0:3]       # X must be estimate of state, not true state

def Ht(X):
    ''' measurement partials '''
    return np.eye(3)    # X must be estimate of state, not true state  

def R(H, t):
    '''
    measurement covariance matrix

    Input:
        - H; dilution of precision matrix
        - t; current time
    '''

    # Clock bias error
    #  GPS SPS spec: UTCOE < 30ns 95% - 15.3ns 1-sigma
    c = 299792458   # m/s, speed of light
    var = (t * c)**2*hvar(t*700)
    # print(sqrt(vclock))

    # Receiver precision
    #  Can't make Tc*d too fast for true receiver; Tc*d = 10.23 MHz
    #  is a good choice as that is the current P(Y) chipping rate
    Tc = 1/(5.115e6)    # Hz, 5.115 Mcps for LunaNet signal 2
    d = 0.1             # correlation spacing [0.1,1] (lower reduces multipath)
    # Time signal is tracked, larger = better but if user is moving quickly
    # might start confusing measurements
    T = 1.00            # s, averaging time
    CN0 = 42            # dB-Hz, signal power over noise power spectral density
    var += (c*Tc)**2 * (d / (4*T*CN0))

    # Node uncertainty
    #  Assuming each node tracks its position autotrueously, it will have some 
    #  associated covariance which is propagated onto the user
    rss = 13.9              # 3-sigma of RSS position errors
    var += (rss / 3)**2     # mean corresponds to the 50th pctl., 0.675std

    # Multipath
    #  Typical multipath ranges from 1m in a benign environment to 5m in a highly
    #  reflective environment for GPS [Misra and Enge (2006), p177]
    var += (1)**2

    return H * var  # measurement variance at given time step

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
    xstar = np.array([0, 0, -r])
    x_true = np.zeros((m,n))
    x_nom  = np.zeros((m,n))
    y = np.zeros((l,n))

    # Dilution of precision
    DOP = np.zeros((l,l,n))
    t = np.linspace(0, 24*60*60, n)    # time steps (in seconds)

    for j in range(n):
        # Compute true trajectory
        x_true[:,j] = phi_true(w, t[j]) @ x0

        # Compute DOP at each step
        visible, _ = findVisibleSats(x_true[:,j], sats, j, elev=10)
        H = np.diag(computeDOP(np.array(x0 / 1000), visible, j))
        DOP[:,:,j] = np.diag(H[0:3])

        # Compute measurement
        y[:,j] = Y(x_true[:,j], R(DOP[:,:,j], t[j]))

    #np.savetxt('test/true.csv', x_true, delimiter=',')

    # run batch filter
    with Batch(t, xstar, x_true, lambda t: phi(A(w), t), y, G, Ht,
               lambda z: R(DOP[:,:,np.where(t == z)[0][0]], z)
               ) as batch:
        fig = plt.figure()
        ax = plt.axes()
        
        iter = 1
        for i in range(iter):
            batch.evaluate()
            mc, stat = batch.plot(ax, last=True if i == iter - 1 else False)
            # update initial guess
            print(batch.x[:,-1])
            #batch.x0 = batch.x[:,-1]
        #np.savetxt('test/est.csv', batch.x, delimiter=',')

        ax.grid()
        # ax.set_ylim(bottom=0, top=65)
        ax.set_xlabel("Samples")
        ax.set_ylabel("Error (m)")
        ax.set_title(f"RMS Position Uncertainty (Batch filter, GPS geometry)")
        ax.legend(handles=[mc, stat])
        plt.show()