#!/usr/bin/python
import os, sys
import h5py
import glob
from sys import argv
from argparse import ArgumentParser
import numpy as np
current_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(current_path, "..","..","tools"))
sys.path.append(lib_path)

#import Edog_runner as edog_runner
import Simulation
import PlottingTools as plt
import matplotlib.pyplot as mplt


parser = ArgumentParser()
parser.add_argument("config_file", default=None)
parser.add_argument("--run_edog", default=False)
parser.add_argument("--data_path", default=None)
args = parser.parse_args()

app_name = "spatialSummation"
data_path = args.data_path

if args.run_edog==True:
    print "Running edog..."
    config_file = args.config_file
    run_id = edog_runner.run_edog(app_name, config_file)
    data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "Data", run_id)
    print data_path


if data_path==None:
    print "Data path not given...."

data_files = glob.glob1(data_path,'*.h5')
num_data_files = len(data_files)
sim = []

for i in range(num_data_files):
    data_file = os.path.join(data_path, data_files[i])
    f = h5py.File(data_file, "r")
    sim.append(Simulation.Simulation(f))


exp = sim[0]

#Plotting
# data = [
#      [exp.stimulus["spatioTemporal"], "Stimulus"]
#     ,[exp.ganglion["response"]["spatioTemporal"], "Ganglion cell response"]
#     ,[exp.ganglion["impulseResponse"]["spatioTemporal"], "Ganglion cell impulse response"]
#     ]
# plt.animateImshowPlots(data,exp.dt, colorbar = True, save_animation = False, animation_name = "rat")
# mplt.show()



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
fig.savefig(os.path.join(data_path, "rat_cellResponse.pdf"))
