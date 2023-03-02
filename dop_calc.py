from Parsing import *
import numpy as np

if __name__ == "__main__":
    sats = parseGmatData("gmat/GPS_20161231_24sat_1day.txt", gmatReport=True)
    n = sats[0].end
    dop = np.zeros((n,))
    r = 6357.       # radius of moon
    nsats = []

    for t in range(0,n):    # find PDOP at each time step
        visible = findVisibleSats(sats, t, elev=5)
        nsats.append(len(visible))
        H = computeDOP(np.array([0.,0.,-r]), visible, t)
        dop[t] = np.sqrt(H[0,0] + H[1,1] + H[2,2])

    t = np.linspace(0,24,n)

    # Plot position dilution of precision over time
    fig = plt.figure()
    ax1 = plt.subplot(211)
    ax1.plot(t, nsats)
    ax1.grid()
    ax1.set_ylabel("Number of visible sats")
    ax1.set_title("Visible Satellites and PDOP for Earth South Pole")
    ax2 = plt.subplot(212)
    ax2.plot(t, dop)
    ax2.set_ylim(bottom=0, top=4)   # cap y axis so plot is readable
    ax2.grid()
    ax2.set_xlabel("Time (hrs)")
    ax2.set_ylabel("PDOP")
    plt.show()