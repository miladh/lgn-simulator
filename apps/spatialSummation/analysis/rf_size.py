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
sims, output_dir = get_simulations.get_simulation_environment(config_file, False)


# Analysis: --------------------------------------------------------------------
exp = sims[0]

cells = ["ganglion"]
idx = exp.num_points * 0.5
idy = idx
t = exp.time_vec

fig = mplt.figure(figsize=(16,12))
n = 0
for c, cell in enumerate(cells):
    responses = []
    n+=1

    ax  = fig.add_subplot(len(cells), 2, n)
    responses.append([exp.singleCellTemporalResponse(cell, idx, idy), "ganglion"])

    ax.plot(t, responses[-1][0], label = responses[-1][1],zorder=3)
    ax.set_xlabel("t[s]", fontsize= 16)
    ax.set_ylabel("Response",fontsize= 16)
    ax.set_title(cell, fontsize = 16)
    ax.legend(frameon=False)

    n+=1
mplt.tight_layout()
mplt.show()
# fig.savefig(os.path.join(output_dir, "rat_cellResponse.png"))
