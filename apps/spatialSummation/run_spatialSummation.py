#!/usr/bin/python
import os, sys
import h5py
import glob
from sys import argv
from argparse import ArgumentParser
current_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(current_path, "..","..","tools"))
sys.path.append(lib_path)

#import Edog_runner as edog_runner
import Simulation
import PlottingTools as plt
import matplotlib.pyplot as mplt


parser = ArgumentParser()
parser.add_argument("run_edog", default=False)
parser.add_argument("config_file", default=None)
parser.add_argument("data_path", default=None)
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
data = [
     [exp.stimulus["spatioTemporal"], "Stimulus"]
    ,[exp.ganglion["response"]["spatioTemporal"], "Ganglion cell response"]
    ,[exp.ganglion["impulseResponse"]["spatioTemporal"], "Ganglion cell impulse response"]
    ]


plt.animateImshowPlots(data,exp.dt, colorbar = True, save_animation = False, animation_name = "rat")
mplt.show()
