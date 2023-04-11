# library imports
import spiceypy as spice
from autograd import grad, jacobian
from scipy.linalg import expm

# local imports
from .Filter import *
from .Dynamics import integrate


class ExtendedKalman(Filter):
    '''
    Runs a posteriori extended Kalman filter. 
    '''
    def __init__(self, t, x0, vx0, x_true, f, u, Y, Ht, Q, R):
        '''
        Input:
         - t; time steps associated with data
         - x0; state estimate -- near true position, not used w/ a priori (m,)
         - vx0; covariance associated with initial state estimate (m,m)
         - x_true; true user trajectory (m,n)
         - f; state function, accepts 3 args (dt, x, u), returns (m,)
         - u; state input, returns (k,)
         - Y; measurements at each time step (l,n)
         - Ht; measurement model partials, accepts 1 arg (state), returns (l,m)
         - Q; state covariance matrix, accepts 1 arg (dt), returns (m,m)
         - R; measurement covariance matrix, accepts 1 arg (time), returns (l,l)
        '''
        super().__init__(t, x0, x_true, Y, Ht, R)   # initialize parent class

        # Track state, residual, estimate covariance, and Kalman gain
        self.f = f
        self.u = u
        self.Q = Q
        self.vx0 = vx0

        self.x[:,0] = x0
        self.P[:,:,0] = vx0
        self.K = np.zeros((self.m,self.l,self.n))


    def __exit__(self, ex_type, ex_val, ex_traceback):
        ''' Catch exit errors '''
        spice.kclear()              # Unload kernels
        # proceed with normal exit
        super().__exit__(ex_type, ex_val, ex_traceback)

    def evaluate(self, verbose=False):
        '''
        Iterate through sat data and determine estimated state at each time step
        '''
        start = timeit.default_timer()
        looptime = []
        self._i = 0 # reset
        self.step() # start at time step 1

        while self.i < self.n:
            loopstart = timeit.default_timer()
            ti = self.t[self.i]
            dt = ti - self.t[self.i-1]

            # Get local variables to reduce evaluations
            x = self.x[:,self.i-1]                              # x_{k-1|k-1}
            fA = jacobian(lambda c: self.f(dt, c, self.u))      # dF_{k}/dt
            A = expm(fA(x) * dt)                                # F_{k}
            # A = np.zeros(np.shape(A))
            z = self.Y[:,self.i]
            P = self.P[:,:,self.i-1]                            # P_{k-1|k-1}
            R = self.R(ti)
            R[R == float('inf')] = sys.float_info.max           # catch 'inf'

            # Predict step
            # Compute the new state \hat{x_{k|k-1}}
            x = integrate(lambda a, b: self.f(a, b, self.u), [0, dt], x)[:,-1]
            P = A @ P @ np.transpose(A) + self.Q(dt)            # \hat{P_{k|k-1}}

            # Update step
            S = self.Ht @ (P @ np.transpose(self.Ht)) + R       # S_{k}
            self.K[:,:,self.i] = P @ (np.transpose(self.Ht) @ np.linalg.inv(S))  # K_{k}
            IKH = np.eye(self.m) - self.K[:,:,self.i] @ self.Ht # (I - K_{k} * H_{k})
            x = IKH @ x + self.K[:,:,self.i] @ z                # x_{k|k}
            p = x[0:3]
            v = x[3:6]
            if np.dot(p, v) < 0: v = v - np.dot(v, p) / np.dot(p, p) * p
            self.x[:,self.i] = np.concatenate((p, v))
            self.P[:,:,self.i] = IKH @ P                        # P_{k|k}
            self.y[:,self.i] = z - self.Ht @ self.x[:,self.i]   # y_{k|k}

            self.step()             # next time step
            looptime.append(timeit.default_timer() - loopstart)

        stop = timeit.default_timer()

        if verbose:
            print("\n########## Runtime Statistics  ##########")
            print("Evaluation time:\t{:.5f}\ts".format(stop-start))
            print("Average loop time:\t{:.5e}\ts".format(sum(looptime) / len(looptime)))
            print("#########################################\n")

