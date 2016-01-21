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
sims, output_dir = get_simulations.get_simulation_environment(config_file, record=True)


# Analysis: --------------------------------------------------------------------
# exp = sims[0]

responses = []
for exp in sims:
    idx = exp.num_points * 0.5
    idy = idx
    t = exp.time_vec
    print np.mean(exp.singleCellTemporalResponse("ganglion", idx, idy))
    responses.append(np.mean(exp.singleCellTemporalResponse("ganglion", idx, idy)))

mplt.plot(responses)
fig = mplt.figure(figsize=(16,12))
# set_xlabel("t[s]", fontsize= 16)
# set_ylabel("Response",fontsize= 16)
# mplt.tight_layout()
mplt.show()
fig.savefig(os.path.join(output_dir, "rat_cellResponse.png"))
