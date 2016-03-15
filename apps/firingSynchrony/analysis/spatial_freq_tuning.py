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
import get_simulations
import plotting_tools as plt


parser = ArgumentParser()
parser.add_argument("sim_ids", help = "simulation ids")
parser.add_argument("record", help = "record results", type = int)
args = parser.parse_args()
sim_ids = args.sim_ids
record = args.record

sims, output_dir=get_simulations.get_simulation_environment(sim_ids, record=record)

# Analysis: --------------------------------------------------------------------
cell_pos_x = 0.5
cell_pos_y = cell_pos_x

data = {"ganglion": {"K": np.zeros(len(sims)),
                    "responses": np.zeros(len(sims)) }
        ,"relay":   {"K": np.zeros(len(sims)),
                            "responses": np.zeros(len(sims))}
        ,"cortical":   {"K": np.zeros(len(sims)),
                                    "responses": np.zeros(len(sims)) }}
fig = mplt.figure(figsize=(8,6))
for cell in data:
    for j, exp in enumerate(sims):
        idx = exp.integrator.nPointsSpatial * cell_pos_x
        idy = exp.integrator.nPointsSpatial * cell_pos_y
        data[cell]["K"][j] = exp.stimulus.spatialFreq
        data[cell]["responses"][j] = np.mean(
        exp.singleCellTemporalResponse(cell, idx, idy))


# Plot:
mplt.plot(data["ganglion"]["K"], data["ganglion"]["responses"], "o--r", label = "ganglion")
mplt.plot(data["relay"]["K"], data["relay"]["responses"], ",-b", label = "relay")
mplt.plot(data["cortical"]["K"], data["cortical"]["responses"], "g--", label = "relay")

# mplt.xlim(0., 1.)
mplt.xlabel("Spatial freq", fontsize= 16)
mplt.ylabel("Response",fontsize= 16)
mplt.tight_layout()
mplt.legend(loc=1)
mplt.show()
fig.savefig(os.path.join(output_dir, "spatial_freq_tuning.png"))
