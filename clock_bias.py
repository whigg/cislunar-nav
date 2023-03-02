import matplotlib.pyplot as plt
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
    return q1*tau**(-1) + 1/6*q2*tau + 11/120*q3*tau**3

c = 299792458               # m/s, speed of light
t = range(1, 14*86400)      # 1 day
rms = [x*sqrt(hvar(x))*1.96 for x in t]
rms2 = [x*sqrt(hvar(x, q1=2.53e-23, q2=4.22e-24, q3=1.00e-38))*1.96 for x in t]
t = [x / 86400 for x in t]

fig = plt.figure()
ax = plt.axes()
ax.semilogy(t, rms, label="RAFS")
ax.semilogy(t, rms2, label="USO")
ax.hlines([1e-9], [t[1]], [t[-1]], ['gray'], ['dotted'], label="1ns")
ax.hlines([30e-9], [t[1]], [t[-1]], ['gray'], ['dashed'], label="30ns")
ax.grid()
ax.legend()
ax.set_xlabel("Time (days)")
ax.set_ylabel("95% clock error (s)")
ax.set_title("Oscillator stability over time")
plt.show()

'''
tr = range(100, 86400, 100) # one day in seconds
t = []                      # convert to hours
rafsv = []                  # s^2, variance of RAFS
usov  = []                  # s^2, variance of USO
for i in tr:
    t.append(i / 3600)
    rafsv.append((c * i)**2 * hvar(i))
    usov.append((c * i)**2 * hvar(i, q1=2.53e-23, q2=4.22e-24, q3=1e-38))

# Plot rms clock error over time
fig = plt.figure()
ax = plt.axes()
ax.semilogy(t, [sqrt(x) for x in rafsv], label="RAFS")
ax.semilogy(t, [sqrt(x) for x in usov] , label="USO")
ax.grid()
ax.set_xlabel("Time (hrs)")
ax.set_ylabel("$c^2 \\tau^2 * \sigma^2(\\tau)$ (m)")
ax.set_title("Position RMSE due to error of GPS clocks")
ax.legend()
plt.show()
'''
