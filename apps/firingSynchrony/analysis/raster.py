#!/usr/bin/python
import os, sys
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as mplt

import h5py
from glob import glob


current_path = os.path.dirname(os.path.realpath(__file__))
lib_path = [os.path.abspath(os.path.join(current_path,"../../../tools/sumatraTracking")),
           os.path.abspath(os.path.join(current_path,"../../../tools/analysis"))]
[sys.path.append(path) for path in lib_path]
import Simulation as sim
import plotting_tools as plt

outputFilePath =  "/home/milad/Dropbox/projects/lgn/code/lgn-simulator/apps/firingSynchrony/firingSynchrony.h5"
outputFile = glob(outputFilePath)[0]
f = h5py.File(outputFile, "r")
exp = sim.Simulation(f)

# Analysis: --------------------------------------------------------------------

data = [
 [exp.stimulus.spatioTemporal, "Stimulus"]
,[exp.ganglion.response["spatioTemporal"], "Ganglion cell response"]
,[exp.ganglion.impulseResponse["spatioTemporal"], "Ganglion cell impulse response"]
,[exp.relay.response["spatioTemporal"], "Relay cell response"]
,[exp.relay.impulseResponse["spatioTemporal"], "Relay cell impulse response"]
,[exp.cortical.response["spatioTemporal"], "Cortical cell response"]
,[exp.cortical.impulseResponse["spatioTemporal"], "Cortical impulse response"]
]

idx = exp.integrator.nPointsSpatial/2
idy = idx

# Raster plots
spike_train_fb = exp.spikeTrain("relay", idx, idy, num_trails = 10)
spike_train_no_fb = exp.spikeTrain("ganglion", idx, idy, num_trails = 10)
spike_train_cor = exp.spikeTrain("cortical", idx, idy, num_trails = 10)


f = mplt.figure(figsize = (8,12))
plt.raster(spike_train_fb, ax = f.add_subplot(3, 1, 1), title="With feedback")
plt.raster(spike_train_no_fb, ax = f.add_subplot(3, 1, 2), title="Without feedback")
plt.raster(spike_train_cor, ax = f.add_subplot(3, 1, 3), title="Cortical")
mplt.tight_layout()


# impulse response temporal plots
mplt.figure()
impresC = exp.temporalImpulseResponse("cortical", idx, idy)
impresR = exp.temporalImpulseResponse("relay", idx, idy)
mplt.plot(exp.integrator.timeVec, impresR, '-or', label="Relay")
mplt.plot( exp.integrator.timeVec, impresC, '-ob', label="cortical")
mplt.legend()

# Animation
plt.animateImshowPlots(data, exp.integrator.temporalResolution, colorbar = True,
save_animation = False, animation_name = "rat")
mplt.show()
