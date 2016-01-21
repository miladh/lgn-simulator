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
import PlottingTools as plt

parser = ArgumentParser()
parser.add_argument("config_file", default=None)
args = parser.parse_args()
config_file = args.config_file


sims, output_dir=get_simulations.get_simulation_environment(config_file, record=False)


# Analysis: --------------------------------------------------------------------
cell_pos_x = np.linspace(0.5,0.7,5)
cell_pos_y = cell_pos_x

responses = np.zeros([len(cell_pos_x), len(sims)])
fig = mplt.figure(figsize=(8,6))
for i, (x, y) in enumerate(zip(cell_pos_x, cell_pos_y)):
    for j, exp in enumerate(sims):
        idx = exp.num_points * x
        idy = exp.num_points * y
        t = exp.time_vec
        responses[i,j] = np.mean(exp.singleCellTemporalResponse("ganglion", idx, idy))

    label = "{0:.2f}".format(cell_pos_x[i]) + "," + "{0:.2f}".format(cell_pos_y[i])
    mplt.plot(responses[i], "o-", label = label)


mplt.xlabel("Spot diameter", fontsize= 16)
mplt.ylabel("Response",fontsize= 16)
mplt.tight_layout()
mplt.legend(loc=4)
mplt.show()
fig.savefig(os.path.join(output_dir, "rat_cellResponse.png"))
