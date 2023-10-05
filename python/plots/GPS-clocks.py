import matplotlib.pyplot as plt

if __name__ == "__main__":
    clock = lambda t, a0, a1, a2: (a0 + a1*t + a2*(t**2)) * 1e9

    time = range(1, 86400*14, 100)   # two weeks
    # prn1 = [clock(t, 0.469126738608e-03,-0.100044417195e-10, 0.) for t in time]
    # prn2 = [clock(t,-0.647393986583e-03,-0.113686837722e-11, 0.) for t in time]
    # prn3 = [clock(t,-0.611455179751e-04,-0.138697942020e-10, 0.) for t in time]
    # prn4 = [clock(t,-0.196907203644e-03, 0.159161572810e-11, 0.) for t in time]
    # prn5 = [clock(t,-0.663353130221e-04,-0.136424205266e-11, 0.) for t in time]
    prn1 = [clock(t, 0.,-0.100044417195e-10, 0.) for t in time]
    prn2 = [clock(t, 0.,-0.113686837722e-11, 0.) for t in time]
    prn3 = [clock(t, 0.,-0.138697942020e-10, 0.) for t in time]
    prn4 = [clock(t, 0., 0.159161572810e-11, 0.) for t in time]
    prn5 = [clock(t, 0.,-0.136424205266e-11, 0.) for t in time]
    time = [t / 86400 for t in time]

    fig = plt.figure()
    ax = plt.axes()
    ax.plot(time, prn1, label="PRN1")
    ax.plot(time, prn2, label="PRN2")
    ax.plot(time, prn3, label="PRN3")
    ax.plot(time, prn4, label="PRN4")
    ax.plot(time, prn5, label="PRN5")
    ax.grid()
    ax.legend()
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Error from UTC (ns)")
    ax.set_title("Broadcast Clock Corrections (PRN 1-5, 01-JAN-2022 00:00:00)")
    plt.show()


