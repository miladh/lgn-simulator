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
import colormaps as cmaps


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
ds = 63

Nt = sims[0].integrator.nPointsTemporal
t_vec = sims[0].integrator.timeVec

data = {"ganglion": {"orientation": np.zeros(len(sims)),
                    "R_a": len(sims)*[[]],
                    "R_b": len(sims)*[[]]}
        ,"relay":   {"orientation": np.zeros(len(sims)),
                    "R_a": len(sims)*[[]],
                    "R_b": len(sims)*[[]]}
        ,"cortical": {"orientation": np.zeros(len(sims)),
                    "R_a": len(sims)*[[]],
                    "R_b": len(sims)*[[]]}
}

corr = {"ganglion": {"orientation": np.zeros(len(sims)),
                    "cross_corr": np.zeros([2*Nt-1,len(sims)])}
        ,"relay":   {"orientation": np.zeros(len(sims)),
                    "cross_corr": np.zeros([2*Nt-1, len(sims)])}
        ,"cortical": {"orientation": np.zeros(len(sims)),
                    "cross_corr": np.zeros([2*Nt-1, len(sims)])}
}


for cell in data:
    for j, exp in enumerate(sims):
        idx = exp.integrator.nPointsSpatial * cell_pos_x
        idy = exp.integrator.nPointsSpatial * cell_pos_y

        data[cell]["orientation"][j] = exp.stimulus.orientation
        data[cell]["R_a"][j].append(exp.singleCellTemporalResponse(cell, idx-ds, idy))
        data[cell]["R_b"][j].append(exp.singleCellTemporalResponse(cell, idx+ds, idy))

        corr[cell]["orientation"][j] = exp.stimulus.orientation
        corr[cell]["cross_corr"][:,j] = np.correlate(\
        exp.singleCellTemporalResponse(cell, idx-ds, idy), \
        exp.singleCellTemporalResponse(cell, idx+ds, idy), "full")

# Plot:
extent = [corr["ganglion"]["orientation"][0], corr["ganglion"]["orientation"][-1],
          -t_vec[-1]/2.,t_vec[-1]/2. ]

vmin = (corr["ganglion"]["cross_corr"]).min()
vmax = (corr["ganglion"]["cross_corr"]).max()
f, ax = mplt.subplots(3, figsize = (8,12), sharex= True)


ax[0].imshow(corr["ganglion"]["cross_corr"], extent= extent,
            vmin=vmin, vmax=vmax,
            cmap=cmaps.viridis, origin="lower", aspect="auto",
            interpolation="nearest")
ax[0].set_ylabel("Cross correlation",fontsize= 16)


ax[1].imshow(corr["relay"]["cross_corr"], extent= extent,
            cmap=cmaps.viridis, origin="lower", aspect="auto",
            vmin=vmin, vmax=vmax,
            interpolation="nearest")
ax[1].set_ylabel("Cross correlation",fontsize= 16)

im = ax[2].imshow(corr["cortical"]["cross_corr"], extent= extent,
            vmin=vmin, vmax=vmax,
            cmap=cmaps.viridis, origin="lower", aspect="auto",
            interpolation="nearest")
ax[2].set_ylabel("Cross correlation",fontsize= 16)
ax[2].set_xlabel(r"Orientation [$\theta$]", fontsize= 16)


f.colorbar(im, ax=ax.ravel().tolist())
# f.tight_layout()

mplt.show()
f.savefig(os.path.join(output_dir, "spatial_freq_tuning.png"))
