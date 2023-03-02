# library imports
import timeit
import sys
import numpy as np
from scipy.linalg import expm

# local imports
from Parsing import *

class Batch:
    '''
    Runs a posteriori batch filtering
    '''
    
    def __init__(self, t, x0, x_true, phi, Y, G, Ht, R):
        '''
        Input:
         - t; time steps associated with data
         - x0; state estimate -- near true position, not used w/ a priori
         - x_true; trueinal user trajectory (m,n)
         - phi; state transition matrix function, accepts 1 arg (time), returns (m,m)
         - Y; measurements at each time step (l,n)
         - G; measurement model, accepts 1 arg (state), returns (l,)
         - Ht; measurement model partials, accepts 1 arg (state), returns (l,m)
         - R; measurement covariance matrix, accepts 1 arg (time), returns (l,l)
        '''

        # dimensions
        self.l = np.shape(Y)[0]         # dimension of measurements
        self.m = np.shape(x_true)[0]     # dimension of state
        self.n = np.shape(x_true)[1]     # number of data points

        # vectors and matrices
        self.y = np.zeros((self.l,self.n))          # position measurements
        self.P = np.zeros((self.m,self.m,self.n))   # covariance matrix
        self.x = np.zeros((self.m,self.n))          # state estimate
        self.x0 = x0                    # nominal initial state
        self.x_true = x_true            # true state
        self.phi = phi                  # state transition matrix function
        self.Y = Y                      # observed measurements          
        self.G = G                      # measurement model function
        self.Ht = Ht                    # measurement partials
        self.R = R                      # measurement covariance matrix function
        self.t = t                      # list of time steps

        self._i = 0                     # initialize sim step

    def __enter__(self):
        ''' Support for with statements '''
        return self

    def __exit__(self, ex_type, ex_val, ex_traceback):
        ''' Catch exit errors '''
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

    def evaluate(self, verbose=False):
        '''
        Iterate through sat data and determine estimated state at each time step
        '''
        start = timeit.default_timer()
        looptime = []
        self._i = 0     # start at zero

        HTWH = np.zeros((self.m, self.m))
        HTWy = np.zeros((self.m,))

        while self.i < self.n:  # start at time step 0
            loopstart = timeit.default_timer()

            ti = self.t[self.i]
            # Do math
            xi = self.phi(ti) @ self.x0
            R = self.R(ti)
            # catch inf or negative vals
            R[(R == float('inf')) | (R < 0)] = sys.float_info.max
            W = np.linalg.inv(R)
            H = self.Ht(xi) @ self.phi(ti)

            self.y[:,self.i] = self.Y[:,self.i] - self.G(xi)
            HTWH = HTWH + H.T @ W @ H
            HTWy = HTWy + H.T @ W @ self.y[:,self.i]

            self.P[:,:,self.i] = np.linalg.inv(HTWH)
            self.x[:,self.i]   = self.x0 + np.linalg.solve(HTWH, HTWy)

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
        t = range(self.n)   # get number of datapoints
        e = [0] * self.n
        std = [0] * self.n
        for j in range(0,self.n):
            diff = self.x[:,j] - self.x_true[:,0]
            var = np.diag(self.P[:,:,j])
            e[j] = np.sqrt(np.dot(diff, diff))              # rms pos error
            std[j] = 3 * np.sqrt(var[0] + var[1] + var[2])  # std of est.

        # Plot RMS position error
        mc, = ax.plot(t, e, color='gray', linestyle='dotted', linewidth=1, label='Monte-Carlo run')
        stat, = ax.plot(t, std, color='green', label='3σ error') if last else (None,)
        return mc, stat

