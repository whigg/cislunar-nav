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
    for file in os.listdir(dir):    # create list of files in directory
        f = os.path.join(dir, file)

        if os.path.isfile(f):
            files.append(f)

    np.random.seed(4020)            # fixed seed

    # plotting
    default_cycler = (cycler(color=['navy','darkgreen','blue','green','cyan','lime']) + cycler(linestyle=['-','--','-','--','-','--']))
    # default_cycler = (cycler(color=['b','g','r','c','m','darkorange']) + cycler(linestyle=['-','--',':','-.','-','--']))
    # default_cycler = (cycler(color=['r','c','m','b','g']) + cycler(linestyle=[':','-.','-','-','--']))
    plt.rc('axes', prop_cycle=default_cycler)
    fig = plt.figure(figsize=(10,2.5))
    ax = plt.axes()
    labels = ['Khon1-4 (100m)', 'Khon1-6 (100m)', 'Khon1-4 (10m)', 'Khon1-6 (10m)', 'Khon1-4 (1m)', 'Khon1-6 (1m)']
    plot = PlotSettings(ax, semilog=True)

    # create and run simulation
    sim = NavSim(files, plot, user='GS', alg='BATCH', debug=True)
    sim.odErr = 100
    sim.run()
    sim.odErr = 10
    sim.run()
    sim.odErr = 1
    sim.run()

    # generate plot
    ax.grid()
    ax.set_ylim(bottom=0, top=40)
    ax.set_xlim(left=0, right=24)   # bound to actual limits
    ax.set_xlabel("Time (hrs)")
    ax.set_ylabel("Error (m)")
    ax.set_title(f"Ground Station 3σ RMS Position Uncertainty, Batch Filter")

    # position legend outside and to right of plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(labels, loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()