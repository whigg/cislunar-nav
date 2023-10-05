from .Parsing import *
from math import sqrt

def hvar(tau, q1=3.7e-24, q2=1.87e-33, q3=7.56e-59):
    '''
    Returns the Hadamard variance of a given oscillator
    rms error after tau is tau*sqrt(hvar(tau))
        tau - time span
        q1  - clock bias (default RAFS)
        q2  - clock rate (default RAFS)
        q3  - clock acceleration (default RAFS)
    '''
    if tau == 0: return 0
    return q1*tau**(-1) + 1/6*q2*tau + 11/120*q3*tau**3

if __name__ == "__main__":
    var = 0

    # Dilution of precision
    #sats = parseGmatData("data/SattPositionsBaseline.txt")
    sats = parseGmatData("data/const_eph/4_MoonOrb.txt", gmatReport=True)
    n = sats[0].end
    r = 1737.4      # radius of moon
    dop = np.zeros((n,))
    pos = np.array([0, 0, -r])

    for t in range(n):
        visible, _ = findVisibleSats(pos, sats, t, elev=0)
        H = computeDOP(pos, visible, t)
        dop[t] = sqrt(abs(H[0,0]) + abs(H[1,1]) + abs(H[2,2]))

    # Clock bias error
    #  GPS SPS spec: UTCOE < 30ns 95% - 15.3ns 1-sigma
    c = 299792458   # m/s, speed of light
    var = np.array([(t*700 * c)**2*hvar(t*700) for t in range(n)])
    # print(sqrt(vclock))

    # Receiver precision
    #  Can't make Tc*d too fast for true receiver; Tc*d = 10.23 MHz
    #  is a good choice as that is the current P(Y) chipping rate
    Tc = 1/(5.115e6)    # Hz, 5.115 Mcps for LunaNet signal 2
    d = 0.1             # correlation spacing [0.1,1] (lower reduces multipath)
    # Time signal is tracked, larger = better but if user is moving quickly
    # might start confusing measurements
    T = 0.05            # s, averaging time
    CN0 = 42            # dB-Hz, signal power over noise power spectral density
    var += (c*Tc)**2 * (d / (4*T*CN0))
    print(f"receiver: {sqrt((c*Tc)**2 * (d / (4*T*CN0)))*1.96}\n")

    # Node uncertainty
    #  Assuming each node tracks its position autonomously, it will have some 
    #  associated covariance which is propagated onto the user
    rss = 13.9              # 3-sigma of RSS position errors
    var += (rss / 3)**2     # mean corresponds to the 50th pctl., 0.675std
    print(f"OD uncertainty: {sqrt((rss / 3)**2)*1.96}\n")

    # Multipath
    #  Typical multipath ranges from 1m in a benign environment to 5m in a highly
    #  reflective environment for GPS [Misra and Enge (2006), p177]
    var += (1)**2

    ure = np.sqrt(var) * 1.96

    fig = plt.figure()
    ax = plt.axes()
    ax.plot(np.linspace(0,24,n), dop)
    ax.grid()
    ax.set_ylim(bottom=1, top=8)  # cap y axis so plot is readable
    ax.set_xlim(left=0, right=24)   # bound to actual limits
    ax.set_xlabel("Time (hrs)")
    ax.set_ylabel("95% Position Error (m)")
    ax.set_title("UNE on lunar south pole, 8 satellites")
    plt.show()

    print(min(dop))
    pass