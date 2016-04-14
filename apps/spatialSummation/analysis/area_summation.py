#!/usr/bin/python
import os, sys
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as mplt

current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../../tools")))
import sumatraTracking.get_simulations as smt

parser = ArgumentParser()
parser.add_argument("sim_ids", help = "simulation ids")
parser.add_argument("record", help = "record results", type = int)
args = parser.parse_args()
sim_ids = args.sim_ids
record = args.record

sims = smt.get_simulations(sim_ids)

output_dir = None
if(record):
    output_dir = smt.get_output_dir(sim_ids)

# Analysis: --------------------------------------------------------------------
cell_pos_x = 0.5
cell_pos_y = cell_pos_x

data = {"ganglion": {"spot_diameter": np.zeros(len(sims)),
                    "responses": np.zeros(len(sims)) }
        ,"relay":   {"spot_diameter": np.zeros(len(sims)),
                            "responses": np.zeros(len(sims))}
        ,"cortical":   {"spot_diameter": np.zeros(len(sims)),
                                    "responses": np.zeros(len(sims)) }}
for cell in data:
    for j, exp in enumerate(sims):
        idx = exp.integrator.nPointsSpatial * cell_pos_x
        idy = exp.integrator.nPointsSpatial * cell_pos_y
        data[cell]["spot_diameter"][j] = exp.stimulus.maskSize
        res = exp.singleCellTemporalResponse(cell, idx, idy)
        res = res[np.where(res  >= 0)]
        data[cell]["responses"][j] = np.mean(res)


# Plot:
fig = mplt.figure(figsize=(8,6))
mplt.plot(data["ganglion"]["spot_diameter"], data["ganglion"]["responses"], "o--r", label = "ganglion")
mplt.plot(data["relay"]["spot_diameter"], data["relay"]["responses"], ">--b", label = "relay")
mplt.plot(data["cortical"]["spot_diameter"], data["cortical"]["responses"], "<--g", label = "cortical")

# mplt.ylim(0., 1.)
mplt.xlabel(r"Spot diameter [deg]", fontsize= 16)
mplt.ylabel("Response",fontsize= 16)
mplt.tight_layout()
mplt.legend(loc=1)
mplt.show()
if record : fig.savefig(os.path.join(output_dir, "area_response.png"))


fig = mplt.figure(figsize=(8,6))
mplt.plot(data["ganglion"]["spot_diameter"], data["ganglion"]["responses"]/ (data["ganglion"]["responses"]).max(), "o--r", label = "ganglion")
mplt.plot(data["relay"]["spot_diameter"], data["relay"]["responses"]/(data["relay"]["responses"]).max(), ">--b", label = "relay")
mplt.plot(data["cortical"]["spot_diameter"], data["cortical"]["responses"]/ (data["cortical"]["responses"]).max(), "<--g", label = "cortical")

mplt.xlabel(r"Spot diameter [deg]", fontsize= 16)
mplt.ylabel("Normalized response",fontsize= 16)
mplt.tight_layout()
mplt.legend(loc=1)
mplt.show()
if record : fig.savefig(os.path.join(output_dir, "area_response_normalized.png"))
