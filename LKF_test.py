# local imports
from filters.LinearKalman import *
from URE import *

def A(dt):
    ''' State matrix at time step dt '''
    A = np.concatenate((np.eye(3), np.eye(3)*dt), axis=1)
    B = np.concatenate((np.zeros((3,3)), np.eye(3)), axis=1)
    return np.concatenate((A, B), axis=0)

def B(dt):
    ''' Acceleration state matrix with time step dt '''
    G = np.array([[0.5*dt**2,0,0],[0,0.5*dt**2,0],[0,0,0.5*dt**2]])
    return np.concatenate((G, np.eye(3)*dt), axis=0)

def Y(X, R):
    ''' true measurement '''
    # X is true state and R is measurement covariance
    return np.random.normal(loc=X[0:3], scale=np.sqrt(np.diag(R)))

def Ht(X):
    ''' measurement partials '''
    # X must be estimate of state, not true state
    return np.concatenate((np.eye(3), np.zeros((3,3))), axis=1)

def Q(B, var):
    ''' State covariance matrix at time step i '''
    return B @ np.transpose(B) * var

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
    T = 0.02            # s, averaging time
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

def nbodydyn(sat, t0, t):
    '''
    Get forces on spacecraft from the planets at time step i
    '''
    merc  = spice.spkezr('MERCURY', t0 + t, 'J2000', 'NONE', 'MOON')[0][0:3]
    venus = spice.spkezr('VENUS', t0 + t, 'J2000', 'NONE', 'MOON')[0][0:3]
    earth = spice.spkezr('EARTH', t0 + t, 'J2000', 'NONE', 'MOON')[0][0:3]
    moon  = spice.spkezr('MOON' , t0 + t, 'J2000', 'NONE', 'MOON')[0][0:3]
    jupi  = spice.spkezr('JUPITER t0TER', t0 + t, 'J2000', 'NONE', 'MOON')[0][0:3]
    saturn= spice.spkezr('SATURN BARYCENTER', t0 + t, 'J2000', 'NONE', 'MOON')[0][0:3]
    uranus= spice.spkezr('URANUS BARYCENTER', t0 + t, 'J2000', 'NONE', 'MOON')[0][0:3]
    nept  = spice.spkezr('NEPTUNE BARYCENTER', t0 + t, 'J2000', 'NONE', 'MOON')[0][0:3]
    sun   = spice.spkezr('SUN'  , t0 + t, 'J2000', 'NONE', 'MOON')[0][0:3]
    sat = sat[0:3]  # take only position vector (if p + v is given)

    uMe = spice.bodvrd("MERCURY", "GM", 1)[1][0]
    uV = spice.bodvrd("VENUS", "GM", 1)[1][0]
    uE = spice.bodvrd("EARTH", "GM", 1)[1][0]
    uM = spice.bodvrd("MOON" , "GM", 1)[1][0]
    uJ = spice.bodvrd("JUPITER", "GM", 1)[1][0]
    uSa = spice.bodvrd("SATURN", "GM", 1)[1][0]
    uU = spice.bodvrd("URANUS", "GM", 1)[1][0]
    uN = spice.bodvrd("NEPTUNE", "GM", 1)[1][0]
    uS = spice.bodvrd("SUN"  , "GM", 1)[1][0]

    GM = [uMe, uV, uE, uM, uJ, uSa, uU, uN, uS]

    rMercSat  = sat - merc
    rVenusSat = sat - venus
    rEarthSat = sat - earth
    rMoonSat  = sat - moon
    rJupiSat  = sat - jupi
    rSaturnSat= sat - saturn
    rUranusSat= sat - uranus
    rNeptSat  = sat - nept
    rSunSat   = sat - sun

    R = [rMercSat, rVenusSat, rEarthSat, rMoonSat, rJupiSat, rSaturnSat, rUranusSat, rNeptSat, rSunSat]

    F = 0
    for j,r in enumerate(R):
        F -= GM[j] / np.linalg.norm(r)**3 * r

    return F

if __name__ == "__main__":
    # Data and dimensions
    sats = parseGmatData("gmat/GPS_20161231_24sat_1day.txt", gmatReport=True)
    n = sats[0].end
    l = 3; m = 6; n = 126    
        
    # Constants
    r = 1737400                             # m, radius of moon
    w = 2*np.pi / (27.3217 * 24*60*60)      # rad/s, rotation rate of moon
    # x0 = np.array([np.sin(np.pi/12)*r, 0., -np.cos(np.pi/12)*r, 0, w*np.sin(np.pi/12)*r, 0])
    x0 = np.array([0., 0., -r, 0., 0., 0.])
    vx0 = np.array([15**2, 15**2, 15**2, 0., 0., 0.])
    xstar = np.array([0, 0, -r])
    x_true = np.zeros((m,n))
    x_nom  = np.zeros((m,n))
    y = np.zeros((l,n))

    # Dilution of precision
    DOP = np.zeros((l,l,n))
    t = np.linspace(0, 24*60*60, n)    # time steps (in seconds)

    for j in range(n):
        # Compute true trajectory
        x_true[:,j] = x0

        # Compute DOP at each step
        visible, _ = findVisibleSats(x_true[0:3,j], sats, j, elev=10)
        H = np.diag(computeDOP(np.array(x0[0:3] / 1000), visible, j))
        DOP[:,:,j] = np.diag(H[0:3])

        # Compute measurement
        y[:,j] = Y(x_true[:,j], R(DOP[:,:,j], t[j]))


    with LinearKalman(t, x0, np.diag(vx0), x_true, A, B, y, Ht(0),
                      lambda z: Q(B(z), (1.5e-3)**2),
                      lambda z: R(DOP[:,:,np.where(t == z)[0][0]], z)
        ) as dyn:

        fig = plt.figure()
        ax = plt.axes()
        for i in range(1):
            dyn.evaluate()
            mc, stat = dyn.plot(ax, last=True if i == 0 else False)

        ax.grid()
        ax.set_ylim(bottom=0, top=100)
        ax.set_xlabel("Time (hrs)")
        ax.set_ylabel("Error (m)")
        ax.set_title(f"RMS Position Uncertainty (Kalman filter, GPS geometry)")
        ax.legend(handles=[mc, stat])
        plt.show()