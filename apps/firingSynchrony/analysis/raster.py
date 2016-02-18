#!/usr/bin/python
import os, sys
from argparse import ArgumentParser
from pylab import*
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
N = float(exp.integrator.nPointsSpatial)
idx = exp.integrator.nPointsSpatial/2
idy = idx

#
# # Raster plots------------------------------------------------------------------
# num_trails = 40
# spike_train_fb = exp.spikeTrain("relay", idx, idy, num_trails = num_trails)
# spike_train_no_fb = exp.spikeTrain("ganglion", idx, idy, num_trails = num_trails)
# spike_train_cor = exp.spikeTrain("cortical", idx, idy, num_trails = 10)
#
# f = mplt.figure(figsize = (8,12))
# plt.raster(spike_train_fb, ax = f.add_subplot(3, 1, 1), title="With feedback")
# plt.raster(spike_train_no_fb, ax = f.add_subplot(3, 1, 2), title="Without feedback")
# plt.raster(spike_train_cor, ax = f.add_subplot(3, 1, 3), title="Cortical")
# mplt.tight_layout()
#
#
# # Response plots------------------------------------------------------------------
# ds = 10
# res_r1 = exp.singleCellTemporalResponse("relay", idx, idy)
# res_r2 = exp.singleCellTemporalResponse("relay", idx+ds, idy+ds)
# res_g1 = exp.singleCellTemporalResponse("ganglion", idx, idy)
# res_g2 = exp.singleCellTemporalResponse("ganglion", idx+ds, idy+ds)
# res_c1 = exp.singleCellTemporalResponse("cortical", idx, idy)
# res_c2 = exp.singleCellTemporalResponse("cortical", idx+ds, idy+ds)
#
# f, ax = mplt.subplots(3, figsize = (8,12))
# ax[0].set_title("With feedback")
# ax[0].plot(res_r1, label= str(idx/N) +"," +str(idy/N) )
# ax[0].plot(res_r2, label= str(idx/N+ds/N) +"," +str(idy/N+ds/N))
# ax[0].legend()
#
# ax[1].set_title("Without feedback")
# ax[1].plot(res_g1, label= str(idx/N) +"," +str(idy/N))
# ax[1].plot(res_g2, label= str(idx/N+ds/N) +"," +str(idy/N+ds/N))
# ax[1].legend()
#
# ax[2].set_title("cortical")
# ax[2].plot(res_c1, label= str(idx/N) +"," +str(idy/N))
# ax[2].plot(res_c2, label= str(idx/N+ds/N) +"," +str(idy/N+ds/N))
# ax[2].legend()
# mplt.tight_layout()
#
# # Phase plots------------------------------------------------------------------
# vector_sum_fb = zeros(num_trails)
# vector_sum_no_fb = zeros(num_trails)
# dt = exp.integrator.temporalResolution
# for i in range(num_trails):
#     vector_sum_fb[i] = sum(cos(2*pi*spike_train_fb[i]/dt)**2 + sin(2*pi*spike_train_fb[i]/dt)**2)
#     vector_sum_no_fb[i] =sum(cos(2*pi*spike_train_no_fb[i]/dt)**2 + sin(2*pi*spike_train_no_fb[i]/dt)**2)
#
# f= mplt.figure(figsize = (8,12))
# num_bins = 10
# hist(vector_sum_fb, num_bins, histtype="step", lw=2, color="black", label="With feedback")
# hist(vector_sum_no_fb, num_bins, histtype="step", lw=2, color="red", label="Without feedback")
# mplt.legend()
#
# # impulse response temporal plots-----------------------------------------------
# mplt.figure()
# impresC = exp.temporalImpulseResponse("cortical", idx, idy)
# impresR = exp.temporalImpulseResponse("relay", idx, idy)
# mplt.plot(exp.integrator.timeVec, impresR, '-or', label="Relay")
# mplt.plot( exp.integrator.timeVec, impresC, '-ob', label="cortical")
# mplt.legend()
#
# mplt.show()
