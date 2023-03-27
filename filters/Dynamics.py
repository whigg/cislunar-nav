# library imports
import autograd.numpy as np
from scipy.linalg import expm
import spiceypy as spice
from autograd import grad, jacobian
from scipy.integrate import solve_ivp

# local imports
# if __name__ == "__main__": from os import chdir; chdir('../')
from URE import *

'''
Linear integration of a (6,) pos / vel state with acceleration input
'''
def linA(dt):
    ''' State matrix at time step dt '''
    A = np.concatenate((np.eye(3), np.eye(3)*dt), axis=1)
    B = np.concatenate((np.zeros((3,3)), np.eye(3)), axis=1)
    return np.concatenate((A, B), axis=0)

def linB(dt):
    ''' Acceleration state matrix with time step dt '''
    G = np.array([[0.5*dt**2,0,0],[0,0.5*dt**2,0],[0,0,0.5*dt**2]])
    return np.concatenate((G, np.eye(3)*dt), axis=0)

def linQ(B, var):
    ''' State covariance matrix at time step i '''
    return B @ np.transpose(B) * var

def linU():
    #return np.random.normal(loc=np.zeros((3,)), scale=1.5e-3)
    return np.zeros((3,))

'''
Position measurements using navigation constellation
'''
def Y(X, R):
    ''' true measurement '''
    # X is true state and R is measurement covariance
    return np.random.normal(loc=X[0:3], scale=np.sqrt(np.diag(R)))

def G(X):
    ''' measurement model '''
    return X[0:3]       # X must be estimate of state, not true state

def Ht_3x6(X):
    ''' measurement partials '''
    # X must be estimate of state, not true state
    return np.concatenate((np.eye(3), np.zeros((3,3))), axis=1)

def Ht_3x3(X):
    ''' measurement partials '''
    # X must be estimate of state, not true state
    return np.eye(3)

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

'''
Stationary ground station dynamics (3,) position, rotating w/ moon
'''

def statA(w):
    ''' state dynamics matrix '''
    A11 = np.array([[0, -w, 0],[w, 0, 0],[0, 0, 0]])
    # A1 = np.concatenate((A11, np.zeros((3,3))), axis=1)
    # A2 = np.concatenate((np.zeros((3,3)), A11), axis=1)
    # return np.concatenate((A1,A2), axis=0)
    return A11

def statB(dt):
    ''' Acceleration state matrix with time step dt '''
    return np.array([[0.5*dt**2,0,0],[0,0.5*dt**2,0],[0,0,0.5*dt**2]])

def statPhi(A, t):
    ''' state transition matrix '''
    return expm(A * t)

def statPhi_true(w, t):
    ''' true state transition matrix '''
    Pd = np.array([[np.cos(w*t), -np.sin(w*t), 0],[np.sin(w*t), np.cos(w*t), 0],[0, 0, 1]])
    # P1 = np.concatenate((Pd,np.zeros((3,3))), axis=1)
    # P2 = np.concatenate((np.zeros((3,3)),Pd), axis=1)
    # return np.concatenate((P1,P2), axis=0)
    return Pd

'''
Satellite dynamics
'''
def nbodydyn(sat, t0, t):
    '''
    Get forces on spacecraft from the planets at time step i
    NOTE: SPICE MUST BE INITIALIZED PRIOR TO USE
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

'''
Lunar surface user, continuous dynamics
'''
def surfDyn(x, u, g, r, W):
    '''
    dynamical equations for lunar surface user
    Input
     - x; state vector at current time [pos; vel], (6,)
     - u; input acceleration (walk on planet)
     - g; planetary acceleration
     - r; hard planetary radius
     - W; rotation rate vector of planet
    '''
    dx = np.zeros(np.shape(x))
    dx[0:3] = x[3:6]

    p = x[0:3]
    dir = p / np.linalg.norm(p)
    vrel = x[3:6] - np.cross(W, p)      # velocity relative to moon-fixed frame

    u += 2*np.cross(W, vrel) + np.cross(W, np.cross(W, p)) - dir*g 
    u += dir * max(np.dot(u, -dir), 0) if np.linalg.norm(p) <= r else 0

    dx[3:6] = u
    return dx

def surfInt(dt, x0, u, g, r, W):
    '''
    Integrates surfDyn(...) given the time step dt; same inputs.
    '''
    func = lambda _, x: surfDyn(x, u, g, r, W)
    sol = solve_ivp(func, (0, dt), x0, t_eval=[0, dt], rtol=1e-6, atol=1e-8)
    if sol.success:
        return sol.y[:,-1]
    else:
        print("Integration Error: " + sol.message)


def tester(dt, x0, *args):
    A = np.array([[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
    func = lambda _, x: A @ x
    sol = solve_ivp(func, [0, dt], x0, t_eval=[0, dt])
    return sol.y[:,-1]


if __name__ == "__main__":

    u = np.zeros((3,))
    r = 1737.4
    w = 2*np.pi / (27.3217 * 24 * 60 * 60)
    W = np.array([0,0,w])
    g = 1.625e-3
    v0 = 0
    x0 = np.array([r*np.sin(15*np.pi/180), 0, -r*np.cos(15*np.pi/180), np.cos(15*np.pi/180)*v0, w*r*np.sin(15*np.pi/180), np.sin(15*np.pi/180)*v0])

    xf = surfInt(3600, x0, u, g, r, W)
    print(x0)
    print(xf)