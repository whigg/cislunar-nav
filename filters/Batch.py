# local imports
from filters.Filter import *

class Batch(Filter):
    '''
    Runs a posteriori batch filtering
    '''
    def __init__(self, t, x0, x_true, phi, Y, G, Ht, R):
        '''
        Input:
         - t; time steps associated with data
         - x0; state estimate -- near true position, not used w/ a priori
         - x_true; true user trajectory (m,n)
         - phi; state transition matrix function, accepts 1 arg (time), returns (m,m)
         - Y; measurements at each time step (l,n)
         - G; measurement model, accepts 1 arg (state), returns (l,)
         - Ht; measurement model partials, accepts 1 arg (state), returns (l,m)
         - R; measurement covariance matrix, accepts 1 arg (time), returns (l,l)
        '''
        super().__init__(t, x0, x_true, Y, Ht, R)

        # vectors and matrices
        self.phi = phi                  # state transition matrix function     
        self.G = G                      # measurement model function


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
            R[R == 0] = 0
            try:
                W = np.linalg.inv(R)
            except np.linalg.LinAlgError:
                W = np.zeros(np.shape(R))
            H = self.Ht(xi) @ self.phi(ti)

            self.y[:,self.i] = self.Y[:,self.i] - self.G(xi)

            HTWH = HTWH + H.T @ W @ H
            HTWy = HTWy + H.T @ W @ self.y[:,self.i]

            try:
                self.P[:,:,self.i] = np.linalg.inv(HTWH)
                self.x[:,self.i]   = self.x0 + np.linalg.solve(HTWH, HTWy)
            except:
                pass


            self.step()             # next time step
            looptime.append(timeit.default_timer() - loopstart)

        stop = timeit.default_timer()

        if verbose:
            print("\n########## Runtime Statistics  ##########")
            print("Evaluation time:\t{:.5f}\ts".format(stop-start))
            print("Average loop time:\t{:.5e}\ts".format(sum(looptime) / len(looptime)))
            print("#########################################\n")

