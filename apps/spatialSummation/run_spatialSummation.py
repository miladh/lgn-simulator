#!/usr/bin/python
import os, sys
import h5py
from glob import glob
from sys import argv
from argparse import ArgumentParser
current_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(current_path, "..","..","tools"))
sys.path.append(lib_path)

import Edog_runner as edog_runner



parser = ArgumentParser()
parser.add_argument("run_edog", default=False)
parser.add_argument("config_file", default=None)
parser.add_argument("data_path", default=None)
args = parser.parse_args()

app_name = "spatialSummation"

if args.run_edog:
    print "Running edog..."
    config_file = args.config_file
    edog_runner.run_edog(app_name, config_file)

if data_path==None:
    print "Cannot find data file..."

data_files = glob.glob1(data_path,'*.h5')
num_data_files = len(data_files)
sim = []
for i in range(num_data_files):
    f = h5py.File(data_files[i], "r")
    sim.append(Simulation(f))



#Plotting
