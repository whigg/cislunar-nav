# library imports
import os
from cycler import cycler

# local imports
from .libs.Batch import *
from .libs.Dynamics import *
from .libs.ExtendedKalman import *
from .libs.LinearKalman import *


class NavSim:
    '''
    Parent class for structuring user navigation simulations, including several
    methods (batch, linear and extended)

    TODO: add support for custom users
    '''
    def __init__(self, files, plot, user='GS', alg='EKF', debug=False, **args):
        '''
        Inputs:
         - files: list of data files to simulate
         - plot:  PlotSettings object
         - user:  which of the built-in users to simulate -- GS, ROVER
         - alg:   which algorithm to use -- EKF (default), LKF, BATCH
         - debug: enable debug outputs -- True or False  
        '''

        algs = ['EKF', 'LKF', 'BATCH']
        users = ['GS', 'ROVER']
        if alg not in algs:
            raise TypeError(f'Invalid alg type. Expected one of: {algs}')
        if user not in users:
            raise TypeError(f'Invalid user type. Expected one of: {users}')


        self.files = files
        self.plot = plot
        self.user = user
        self.debug = debug
        self.l = 3      # measurement vector dimension
        self.w = 2*np.pi / (27.3217 * 24*60*60)     # rad/s, mean motion of moon

        if alg == 'BATCH':
            if user == 'ROVER':
                raise ValueError('Incompatible user and algorithm: ROVER and BATCH')
            self.run = self.batchSim
        elif alg == 'EKF':
            self.run = self.EKFSim
        elif alg == 'LKF':
            pass

        self.odErr = 10

    def EKFSim(self):

        rad = 1737.4    # km, radius of moon
        m = 6           # state vector dimension
        ang = 15
        W = np.array([0,0,self.w])
        g = 1.625e-3    # km/s^2, gravity of moon

        for i, file in enumerate(self.files):
            # load data
            sats = parseGmatData(file, gmatReport=True)
            n = sats[0].end

            # Compute true trajectory
            if self.user == 'GS':
                x0 = np.array([rad*np.sin(ang*np.pi/180), 0, -rad*np.cos(ang*np.pi/180), 
                    0., self.w*rad*np.sin(ang*np.pi/180), 0.])
                vx0 = np.array([10**2, 10**2, 10**2, 0.01**2, 0.01**2, 0.01**2])

                t = np.linspace(0, 24*60*60, n)    # time steps (in seconds)
                u = lambda _: np.zeros((3,))
                func = lambda t, x, u: surfDyn(t, x, u, g, rad, W)
                x_true = integrate(lambda t, x: func(t, x, u), t, x0)

                # state covariance
                Q = lambda z: extQ(z, (1e-6)**2)

            elif self.user == 'ROVER':
                v0 = 0.000033
                std_accel = 1e-10    # root Allan variance
                x0 = np.array([rad*np.sin(ang*np.pi/180), 0, -rad*np.cos(ang*np.pi/180), 
                    np.cos(ang*np.pi/180)*v0, self.w*rad*np.sin(ang*np.pi/180), np.sin(ang*np.pi/180)*v0])
                vx0 = np.array([1**2, 1**2, 1**2, 0**2, 0**2, 0**2])

                t = np.linspace(0, 24*60*60, n)    # time steps (in seconds)
                func = lambda t, x, u: surfDyn(t, x, u, g, rad, W)
                randWalk, u = getAccelFuncs(t, 10, 1e-7, std_accel)  
                x_true = integrate(lambda t, x: func(t, x, randWalk), t, x0)

                # state covariance
                Q = lambda z: linQ(statB_6x3(z), (std_accel * 3)**2)

            xstar = np.random.normal(loc=x0, scale=np.sqrt(vx0))
            xstar[0:3] = xstar[0:3] / np.linalg.norm(xstar[0:3]) * rad  # place on surface
            y = np.zeros((self.l,n))

            # Dilution of precision
            DOP = np.zeros((self.l,self.l,n))

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

            # np.savetxt('matlab/true.csv', x_true, delimiter=',')

            # run extended kalman filter
            with ExtendedKalman(t, xstar, np.diag(vx0), x_true, func, u, y, Ht_3x6(0), Q,
                    lambda z: R(DOP[:,:,np.where(t == z)[0][0]], z, rss = self.odErr, receiver=False)/1e6) as dyn:
                
                dyn.evaluate()
                dyn.units = 1e3         # unit conversion for plotting
                dyn.plot(self.plot.ax, mc=self.plot.mc, std=self.plot.std, batch=False, semilog=self.plot.semilog)

                if self.debug:
                    # update initial guess
                    print("Final x: " + str(dyn.x[:,-1]))
                    # print best 3-sigma error
                    var = np.diag(dyn.P[:,:,j])
                    print("3s error: " + str(3 * np.sqrt(var[0] + var[1] + var[2])))

    def batchSim(self):

        rad = 1737400   # m, radius of moon
        m = 3           # state vector dimension

        for i, file in enumerate(self.files):
            # load data
            sats = parseGmatData(file, gmatReport=True)
            n = sats[0].end

            x0 = np.array([rad*np.sin(np.pi/12), 0, -rad*np.cos(np.pi/12)])
            xstar = np.array([0, 0, 0])
            x_true = np.zeros((m,n))
            x_nom  = np.zeros((m,n))
            y = np.zeros((self.l,n))

            # Dilution of precision
            DOP = np.zeros((self.l,self.l,n))
            t = np.linspace(0, 24*60*60, n)    # time steps (in seconds)

            for j in range(n):
                # Compute true trajectory
                x_true[:,j] = statPhi_true(self.w, t[j]) @ x0

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
            with Batch(t, xstar, x_true, lambda t: statPhi(statA(self.w), t), y, G, Ht_3x3,
                    lambda z: R(DOP[:,:,np.where(t == z)[0][0]], z, rss = self.odErr, receiver=False)) as batch:
                
                batch.evaluate()
                batch.plot(self.plot.ax, mc=self.plot.mc, std=self.plot.std, semilog=self.plot.semilog)
                #np.savetxt('test/est.csv', batch.x, delimiter=',')

                if self.debug:
                    # update initial guess
                    print("Final x: " + str(batch.x[:,-1]))
                    # print best 3-sigma error
                    var = np.diag(batch.P[:,:,j])
                    print("3s error: " + str(3 * np.sqrt(var[0] + var[1] + var[2])))

class PlotSettings:
    '''
    wrapper to contain plot settings

    Input:
     - ax:      plot axes, generated by matplotlib.pyplot.axes()
     - mc:      boolean, include Monte-Carlo results in plot?
     - std:     boolean, include 3-sigma standard deviation line in plot?
     - semilog: boolean, make y-axis log scale?
    '''

    def __init__(self, ax, mc=False, std=True, semilog=False):
        self.ax = ax
        self.mc = mc
        self.std = std
        self.semilog = semilog