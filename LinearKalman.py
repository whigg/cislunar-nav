from Parsing import *
from data.Training import *
from URE import *

import spiceypy as spice
import sys
import timeit

class LinearKalman:
    '''
    Runs a posteriori Kalman filter on GMAT simulation data. Used to analyze
    achievable navigation solution for given lunar relay. 

    Input:
     data - string path to GMAT data file for true position and time
     vaccl - float variance of environmental uncertainty
     vmeas - float variance of nav measurements
     x0    - vector of initial position and velocty (6,)
     vx0   - matrix of uncertainty for x0 (6,6)
    '''
    
    def __init__(self, data, vaccl, x0, vx0, gmat=False):
        '''
        Get user data from GMAT simulation in format
         user.X  user.Y  user.Z
         x(0)    y(0)    z(0)
         ...     ...     ...
         x(n)    y(n)    z(n)

        and initialize filter
        '''
        sats = parseGmatData(data, gmatReport=gmat)
        n = sats[0].end
        self.n = n
        self.x0 = x0                # position the user will be at inertially at all times
        self.vx0 = vx0

        # Dilution of precision
        self.DOP = np.zeros((3,3,n))

        for t in range(n):
            visible, _ = findVisibleSats(x0[0:3], sats, t, elev=0)
            H = np.diag(computeDOP(np.array(self.x0[0:3]), visible, t))
            self.DOP[:,:,t] = np.diag(H[0:3])

        self.t = np.linspace(0, 24*60*60, n)    # time steps (in seconds)
        self._i = 0                             # initialize sim step

        self.vaccl = vaccl
        self._i = 0

        # Track state, residual, estimate covariance, and Kalman gain
        self.x = np.zeros((6,self.n))
        self.y = np.zeros((3,self.n))
        self.P = np.zeros((6,6,self.n))
        self.K = np.zeros((6,3,self.n))
        self.e = self.x

        self.P[:,:,0] = vx0

        self.H = np.concatenate((np.eye(3), np.zeros((3,3))), axis=1)


    def __enter__(self):
        ''' Support for with statements '''
        return self

    def __exit__(self, ex_type, ex_val, ex_traceback):
        ''' Catch exit errors '''
        spice.kclear()              # Unload kernels

        if ex_type is not None:
            print("\nExecution type:", ex_type)
            print("\nExecution value:", ex_val)
            print("\nTraceback:", ex_traceback)

    @property
    def i(self):
        return self._i

    def step(self):
        '''
        Increment time step of simulation
        '''
        self._i += 1

    @property
    def F(self):
        '''
        State matrix at time step i
        '''
        t = self.t[self.i] - self.t[self.i-1]
        A = np.concatenate((np.eye(3), np.eye(3)*t), axis=1)
        B = np.concatenate((np.zeros((3,3)), np.eye(3)), axis=1)
        return np.concatenate((A, B), axis=0)

    @property
    def G(self):
        ''' Acceleration state matrix at time step i '''
        t = self.t[self.i] - self.t[self.i-1]
        G = np.array([[0.5*t**2,0,0],[0,0.5*t**2,0],[0,0,0.5*t**2]])
        return np.concatenate((G, np.eye(3)*t), axis=0)

    @property
    def Q(self):
        '''
        State covariance matrix at time step i
        '''
        G = self.G
        return G @ np.transpose(G) * self.vaccl
        # return self.vaccl

    @property
    def R(self):
        H = self.DOP[:,:,self.i]  # current DOP matrix
        t = self.t[self.i]      # current time step
    
        # Clock bias error
        #  GPS SPS spec: UTCOE < 30ns 95% - 15.3ns 1-sigma
        c = 299792458   # m/s, speed of light
        var = (t * c)**2*hvar(t*700)
        # print(sqrt(vclock))

        # Receiver precision
        #  Can't make Tc*d too fast for nominal receiver; Tc*d = 10.23 MHz
        #  is a good choice as that is the current P(Y) chipping rate
        Tc = 1/(5.115e6)    # Hz, 5.115 Mcps for LunaNet signal 2
        d = 0.1             # correlation spacing [0.1,1] (lower reduces multipath)
        # Time signal is tracked, larger = better but if user is moving quickly
        # might start confusing measurements
        T = 0.02              # s, averaging time
        CN0 = 42              # dB-Hz, signal power over noise power spectral density
        var += (c*Tc)**2 * (d / (4*T*CN0))

        # Node uncertainty
        #  Assuming each node tracks its position autonomously, it will have some 
        #  associated covariance which is propagated onto the user
        rss = 13.9              # 3-sigma of RSS position errors
        var += (rss / 3)**2     # mean corresponds to the 50th pctl., 0.675std

        # Multipath
        #  Typical multipath ranges from 1m in a benign environment to 5m in a highly
        #  reflective environment for GPS [Misra and Enge (2006), p177]
        var += (1)**2

        return H * var / 1e6    # measurement variance at given time step (km^2)

    def nbodydyn(self, sat):
        '''
        Get forces on spacecraft from the planets at time step i
        '''
        t = self.t[self.i]
        merc  = spice.spkezr('MERCURY', self.et + t, 'J2000', 'NONE', 'MOON')[0][0:3]
        venus = spice.spkezr('VENUS', self.et + t, 'J2000', 'NONE', 'MOON')[0][0:3]
        earth = spice.spkezr('EARTH', self.et + t, 'J2000', 'NONE', 'MOON')[0][0:3]
        moon  = spice.spkezr('MOON' , self.et + t, 'J2000', 'NONE', 'MOON')[0][0:3]
        jupi  = spice.spkezr('JUPITER BARYCENTER', self.et + t, 'J2000', 'NONE', 'MOON')[0][0:3]
        saturn= spice.spkezr('SATURN BARYCENTER', self.et + t, 'J2000', 'NONE', 'MOON')[0][0:3]
        uranus= spice.spkezr('URANUS BARYCENTER', self.et + t, 'J2000', 'NONE', 'MOON')[0][0:3]
        nept  = spice.spkezr('NEPTUNE BARYCENTER', self.et + t, 'J2000', 'NONE', 'MOON')[0][0:3]
        sun   = spice.spkezr('SUN'  , self.et + t, 'J2000', 'NONE', 'MOON')[0][0:3]
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

        norm = lambda x: np.linalg.norm(x)
        F = 0
        for j,r in enumerate(R):
            F -= GM[j] / norm(r)**3 * r

        return F
        # return self.f.vec(self.i)

    def predict(self, sat):
        '''
        Use the MLPRegressor model to predict next state
        '''
        model = self.mlTools[0]
        xscaler = self.mlTools[1]
        yscaler = self.mlTools[2]

        x_test = np.concatenate((sat[0:3], np.array([self.t[self.i]]))).reshape(-1,1).T
        y_test = yscaler.inverse_transform(model.predict(xscaler.transform(x_test)))[0]
        return y_test

    def evaluate(self, verbose=False):
        '''
        Iterate through sat data and determine estimated state at each time step
        '''
        start = timeit.default_timer()
        looptime = []
        self._i = 0 # reset
        self.x[:,0] = np.random.normal(loc=x0, scale=np.sqrt(np.diag(self.vx0)))
        self.step() # start at time step 1

        while self.i < self.n:
            loopstart = timeit.default_timer()

            # Get local variables to reduce evaluations
            x = self.x[:,self.i-1]      # x_{k-1|k-1}
            F = self.F                  # F_{k}
            P = self.P[:,:,self.i-1]    # P_{k-1|k-1}

            R = self.R
            # catch inf or negative vals
            R[(R == float('inf')) | (R < 0)] = sys.float_info.max
            # Introduce measurement noise with std sqrt(vmeas)
            z = np.random.normal(loc=self.x0[0:3], scale=np.sqrt(np.diag(R)))

            # Predict step
            accel = np.random.normal(loc=np.zeros((3,)), scale=np.sqrt(self.vaccl))

            x = F @ x + self.G @ accel                  # \hat{x_{k|k-1}}
            self.e[:,self.i] = x
            P = F @ P @ np.transpose(F) + self.Q        # \hat{P_{k|k-1}}

            # Update step
            S = self.H @ (P @ np.transpose(self.H)) + R    # S_{k}
            self.K[:,:,self.i] = P @ (np.transpose(self.H) @ np.linalg.inv(S))  # K_{k}
            IKH = np.eye(6) - self.K[:,:,self.i] @ self.H       # (I - K_{k} * H_{k})
            self.x[:,self.i] = IKH @ x + self.K[:,:,self.i] @ z # x_{k|k}
            self.P[:,:,self.i] = IKH @ P                        # P_{k|k}
            self.y[:,self.i] = z - self.H @ self.x[:,self.i]    # y_{k|k}

            self.step()             # next time step
            looptime.append(timeit.default_timer() - loopstart)

        stop = timeit.default_timer()

        if verbose:
            print("\n########## Runtime Statistics  ##########")
            print("Evaluation time:\t{:.5f}\ts".format(stop-start))
            print("Average loop time:\t{:.5e}\ts".format(sum(looptime) / len(looptime)))
            print("#########################################\n")

    def plot(self, ax, last=False):
        # Compute RMS position error
        t = [x / 3600 for x in self.t]   # get time in hours
        e = [0] * self.n
        std = [0] * self.n
        for j in range(0,self.n):
            diff = self.x[0:3,j] - self.x0[0:3]
            var = np.diag(self.P[:,:,j])
            e[j] = np.sqrt(np.dot(diff, diff)) * 1000           # rms pos error
            std[j] = 3 * np.sqrt(var[0] + var[1] + var[2]) * 1000   # std of est.

        # Plot RMS position error
        mc, = ax.plot(t, e, color='gray', linestyle='dotted', linewidth=1, label='Monte-Carlo run')
        stat, = ax.plot(t, std, color='green', label='3σ error') if last else (None,)
        return mc, stat



if __name__ == "__main__":
    # initial conditions
    r = 1737.4                  # radius of moon
    x0 = np.array([0., 0., -r, 0., 0., 0.])
    vx0 = np.array([0.015**2, 0.015**2, 0.015**2, 0., 0., 0.])

    #file = "data/SattPositionsBaseline.txt"
    file = "gmat/GPS_20161231_24sat_1day.txt"
    typ = True

    with LinearKalman(file, (1.5e-6)**2, x0, np.diag(vx0), gmat=typ) as dyn:

        fig = plt.figure()
        ax = plt.axes()
        for i in range(10):
            dyn.evaluate()
            mc, stat = dyn.plot(ax, last=True if i == 9 else False)

        ax.grid()
        ax.set_ylim(bottom=0, top=100)
        ax.set_xlabel("Time (hrs)")
        ax.set_ylabel("Error (m)")
        ax.set_title(f"RMS Position Uncertainty (Kalman filter, GPS geometry)")
        ax.legend(handles=[mc, stat])
        plt.show()
