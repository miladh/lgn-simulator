#!/usr/bin/python
import os, sys
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as mplt

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


# Analysis: --------------------------------------------------------------------
cell_pos_x = 0.5
cell_pos_y = cell_pos_x

types = ["ganglion", "relay"]
responses = np.zeros([len(types), len(sims)])
spot_diameter = np.zeros(len(sims))

fig = mplt.figure(figsize=(8,6))
for i, t in enumerate(types):
    for j, exp in enumerate(sims):
        idx = exp.num_points * cell_pos_x
        idy = exp.num_points * cell_pos_y
        spot_diameter[j] = exp.stimulus.maskSize
        responses[i,j] = np.mean(exp.singleCellTemporalResponse(t, idx, idy))

    mplt.plot(spot_diameter, responses[i], "o-", label = t)

mplt.xlim(0., 1.)
mplt.xlabel("Spot diameter", fontsize= 16)
mplt.ylabel("Response",fontsize= 16)
mplt.tight_layout()
mplt.legend(loc=1)
mplt.show()
fig.savefig(os.path.join(output_dir, "nonlaggedXCells.png"))
