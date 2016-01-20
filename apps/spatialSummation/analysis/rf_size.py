#!/usr/bin/python
import os, sys
from sys import argv
from argparse import ArgumentParser
import h5py
import numpy as np

current_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(current_path,"..", "..","..","tools"))
sys.path.append(lib_path)
import Simulation
import get_simulations
import PlottingTools as plt
import matplotlib.pyplot as mplt


parser = ArgumentParser()
parser.add_argument("config_file", default=None)
args = parser.parse_args()

config_file = args.config_file
# Read config file:-------------------------------------------------------------
sim = get_simulations.get_simulations(config_file)

# Plotting: --------------------------------------------------------------------
exp = sims[0]

nPoints = 5
cells = ["ganglion"]
idx = np.sort(np.random.randint(0, exp.num_points, nPoints))
idy = np.sort(np.random.randint(0, exp.num_points, nPoints))
t = exp.time_vec

fig = mplt.figure(figsize=(16,12))
n = 0
for c, cell in enumerate(cells):
    spikeTrains = []
    responses = []
    n+=1

    ax  = fig.add_subplot(len(cells), 2, n)
    for i,j in zip(idx, idy):
        label = "(" + str('%.1f' % (i/float(exp.num_points,))) + ", " +  str('%.1f' % (j/float(exp.num_points,)))+ ")"
        spikeTrains.append([exp.spikeTrain(cell, i, j, num_trails = 1), label])
        responses.append([exp.singleCellTemporalResponse(cell, i, j), label])

        ax.plot(t, responses[-1][0], label = responses[-1][1],zorder=3)
        ax.set_xlabel("t[s]", fontsize= 16)
        ax.set_ylabel("Response",fontsize= 16)
        ax.set_title(cell, fontsize = 16)
        ax.legend(frameon=False)

    n+=1
    ax = fig.add_subplot(len(cells), 2, n)
    ax = plt.raster(spikeTrains, figsize = (6,6), ax = ax)
    ax.set_title(cell, fontsize = 16)

mplt.tight_layout()
mplt.show()
fig.savefig(os.path.join(output_dir, "rat_cellResponse.png"))
