# library imports
import timeit
import sys
import numpy as np
import spiceypy as spice

# local imports
from Parsing import *

class LinearKalman:
    '''
    Runs a posteriori Kalman filter. 
    '''
    def __init__(self, t, x0, vx0, x_true, A, B, Y, H, Q, R):
        '''
        Input:
         - t; time steps associated with data
         - x0; state estimate -- near true position, not used w/ a priori (m,)
         - vx0; covariance associated with initial state estimate (m,m)
         - x_true; true user trajectory (m,n)
         - A; state dynamics matrix, accepts 1 arg (dt), returns (m,m)
         - B; input dynamics matrix, accepts 1 arg (dt), returns (m,l)
         - Y; measurements at each time step (l,n)
         - H; measurement model partials, accepts 1 arg (state), returns (l,m)
         - Q; state covariance matrix, accepts 1 arg (dt), returns (m,m)
         - R; measurement covariance matrix, accepts 1 arg (time), returns (l,l)
        '''

        # dimensions
        self.l = np.shape(Y)[0]         # dimension of measurements
        self.m = np.shape(x_true)[0]    # dimension of state
        self.n = np.shape(x_true)[1]    # number of data points

        self.x0 = x0                # position the user will be at inertially at all times
        self.vx0 = vx0
        self.t = t                  # time steps (in seconds)

        # Track state, residual, estimate covariance, and Kalman gain
        self.x = np.zeros((self.m,self.n))
        self.A = A
        self.B = B
        self.Y = Y
        self.H = H
        self.Q = Q
        self.R = R
        self.y = np.zeros((self.l,self.n))
        self.P = np.zeros((self.m,self.m,self.n))
        self.P[:,:,0] = vx0
        self.K = np.zeros((self.m,self.l,self.n))
        self.e = self.x

        self._i = 0                 # initialize sim step

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

    def evaluate(self, verbose=False):
        '''
        Iterate through sat data and determine estimated state at each time step
        '''
        start = timeit.default_timer()
        looptime = []
        self._i = 0 # reset
        self.x[:,0] = np.random.normal(loc=self.x0, scale=np.sqrt(np.diag(self.vx0)))
        self.step() # start at time step 1

        while self.i < self.n:
            loopstart = timeit.default_timer()
            ti = self.t[self.i]
            dt = ti - self.t[self.i-1]

            # Get local variables to reduce evaluations
            x = self.x[:,self.i-1]      # x_{k-1|k-1}
            A = self.A(dt)              # F_{k}
            B = self.B(dt)
            z = self.Y[:,self.i]
            P = self.P[:,:,self.i-1]    # P_{k-1|k-1}
            R = self.R(ti)
            # catch inf or negative vals
            R[(R == float('inf')) | (R < 0)] = sys.float_info.max

            # Predict step
            accel = np.random.normal(loc=np.zeros((3,)), scale=1.5e-3)

            x = A @ x + B @ accel       # \hat{x_{k|k-1}}
            self.e[:,self.i] = x
            P = A @ P @ np.transpose(A) + self.Q(dt)            # \hat{P_{k|k-1}}

            # Update step
            S = self.H @ (P @ np.transpose(self.H)) + R         # S_{k}
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
            e[j] = np.sqrt(np.dot(diff, diff))              # rms pos error
            std[j] = 3 * np.sqrt(var[0] + var[1] + var[2])  # std of est.

        # Plot RMS position error
        mc, = ax.plot(t, e, color='gray', linestyle='dotted', linewidth=1, label='Monte-Carlo run')
        stat, = ax.plot(t, std, color='green', label='3σ error') if last else (None,)
        return mc, stat

