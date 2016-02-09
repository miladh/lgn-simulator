#!/usr/bin/python
import os, sys
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as mplt

current_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(current_path,"../../../tools"))
sys.path.append(lib_path)
import Simulation
import get_simulations
import plotting_tools as plt

parser = ArgumentParser()
parser.add_argument("sim_ids", help = "simulation ids")
parser.add_argument("record", help = "record results")
args = parser.parse_args()
sim_ids = args.sim_ids
record = args.record

sims, output_dir=get_simulations.get_simulation_environment(sim_ids, record=record)


# Analysis: --------------------------------------------------------------------
cell_pos_x = np.linspace(0.5,0.7,5)
cell_pos_y = cell_pos_x

responses = np.zeros([len(cell_pos_x), len(sims)])
spot_diameter = np.zeros(len(sims))

fig = mplt.figure(figsize=(8,6))
for i, (x, y) in enumerate(zip(cell_pos_x, cell_pos_y)):
    for j, exp in enumerate(sims):
        idx = exp.num_points * x
        idy = exp.num_points * y
        spot_diameter[j] = exp.stimulus.maskSize
        responses[i,j] = np.mean(exp.singleCellTemporalResponse("ganglion", idx, idy))

    label = "{0:.2f}".format(cell_pos_x[i]) + "," + "{0:.2f}".format(cell_pos_y[i])
    mplt.plot(spot_diameter, responses[i], "o-", label = label)


mplt.xlabel("Spot diameter", fontsize= 16)
mplt.ylabel("Response",fontsize= 16)
mplt.tight_layout()
mplt.legend(loc=4)
mplt.show()
fig.savefig(os.path.join(output_dir, "rat_cellResponse.png"))
