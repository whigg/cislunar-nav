# library imports
import os
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler

# local imports
from filters.NavSim import *

if __name__ == "__main__":
    # Data
    dir = 'data/IM/'                # relative path to data files
    files = []
    # for file in os.listdir(dir):    # create list of files in directory
    #     f = os.path.join(dir, file)

    #     if os.path.isfile(f):
    #         files.append(f)

    files.append(os.path.join(dir, 'Khon2-7_MOON_ME.txt'))
    title = 'Khon 2-7'

    # plotting
    # default_cycler = (cycler(color=['navy','darkgreen','blue','green','cyan','lime']) + cycler(linestyle=['-','--','-','--','-','--']))
    default_cycler = (cycler(color=['b','g','r','c','m','darkorange']) + cycler(linestyle=['-','--',':','-.','-','--']))
    # default_cycler = (cycler(color=['r','c','m','b','g']) + cycler(linestyle=[':','-.','-','-','--']))
    plt.rc('axes', prop_cycle=default_cycler)
    fig, axs = plt.subplots(2, figsize=(10,5))
    labels = ['Monte-Carlo run', '3σ bound']
    plot = PlotSettings(axs[1], semilog=False, mc=True)

    # create and run simulation(s)
    sim = NavSim(files, plot, user='ROVER', alg='EKF', debug=True)
    cases = [38.2]
    for case in cases:
        sim.odErr = case
        # np.random.seed(69)              # fixed seed
        sim.run()

    # generate plot
    axs[1].grid()
    axs[1].set_ylim(bottom=0, top=50)
    axs[1].set_xlim(left=0, right=24)   # bound to actual limits
    axs[1].set_xlabel("Time (hrs)")
    axs[1].set_ylabel("Error (m)")
    axs[1].set_title(f"Rover 3σ RMS Position Uncertainty, EKF")
    axs[1].legend(labels, loc='upper right')

    # # position legend outside and to right of plot
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # ax.legend(labels, loc='center left', bbox_to_anchor=(1, 0.5))

    axs[0].plot(sim.t / 3600, sim.nvis)
    axs[0].grid()
    axs[0].set_ylim(bottom=0, top=5)
    axs[0].set_xlim(left=0, right=24)   # bound to actual limits
    axs[0].set_ylabel('# of Visible Satellites')
    axs[0].set_title('Visible Satellites to User')

    plt.suptitle('24-hour Navigation Simulation (6 Feb 2027 00:00:00)', weight='bold')
    plt.tight_layout()
    plt.show()