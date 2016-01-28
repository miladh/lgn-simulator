#!/usr/bin/python
import os, sys
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as mplt
from scipy import io

current_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(current_path,"..", "..","..","tools"))
sys.path.append(lib_path)
import Simulation
import get_simulations
import plotting_tools as plt

parser = ArgumentParser()
parser.add_argument("config_file", default=None)
args = parser.parse_args()
config_file = args.config_file


sims, output_dir=get_simulations.get_simulation_environment(config_file, record=False)


# Get experimental data: -------------------------------------------------------
exp_data_path = os.path.abspath(os.path.join(current_path,"../../../.." ,"expDATA"))
exp_data_filename = "nonLaggedXCells.mat"
exp_data = io.loadmat(os.path.join(exp_data_path, exp_data_filename ))
relay_exp = [np.array(exp_data["xrelay_data"]).transpose()[0]*2.0,
            np.array(exp_data["yrelay_data"]).transpose()[0]]
ganglion_exp = [np.array(exp_data["xgang_data"]).transpose()[0]*2.0,
               np.array(exp_data["ygang_data"]).transpose()[0]]


# Analysis: --------------------------------------------------------------------
cell_pos_x = 0.5
cell_pos_y = cell_pos_x

data = {"ganglion": {"spot_diameter": np.zeros(len(sims)),
                    "responses": np.zeros(len(sims)) }
        ,"relay":   {"spot_diameter": np.zeros(len(sims)),
                            "responses": np.zeros(len(sims)) }}

for cell in data:
    for j, exp in enumerate(sims):
        idx = exp.num_points * cell_pos_x
        idy = exp.num_points * cell_pos_y
        data[cell]["spot_diameter"][j] = exp.stimulus.maskSize
        data[cell]["responses"][j] = np.mean(
        exp.singleCellTemporalResponse(cell, idx, idy))


# Rescale
data["ganglion"]["spot_diameter"] *= 13.5
data["ganglion"]["responses"] *= 128
data["ganglion"]["responses"] += 36

data["relay"]["spot_diameter"] *= 13.5
data["relay"]["responses"] *= 110
data["relay"]["responses"] += 6


# Plot:
mplt.plot(data["ganglion"]["spot_diameter"], data["ganglion"]["responses"], "r-", label = "ganglion")
mplt.plot(data["relay"]["spot_diameter"], data["relay"]["responses"], "b-", label = "relay")

# Plot experimental data:
mplt.plot(ganglion_exp[0], ganglion_exp[1], "ro", markersize=8)
mplt.plot(relay_exp[0], relay_exp[1], "bo", markersize=8)



# mplt.xlim(0., 1.)
mplt.xlabel(r"Spot diameter [deg]", fontsize= 16)
mplt.ylabel("Response [spikes/s]",fontsize= 16)
mplt.tight_layout()
mplt.legend(loc=1)
mplt.show()
fig.savefig(os.path.join(output_dir, "nonlaggedXCells.png"))
