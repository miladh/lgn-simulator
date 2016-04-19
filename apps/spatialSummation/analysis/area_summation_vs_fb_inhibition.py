#!/usr/bin/python
import os, sys
from argparse import ArgumentParser
import numpy as np
from pylab import*
import operator

current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../../tools")))
import sumatraTracking.io_manager as smt
import analysis.colormaps as cmaps
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
data = {}
for sim in sims:
    d = getattr(getattr(sim, "stimulus"), "maskSize")
    data.setdefault(d, []).append(sim)


def area_summation(sims, cell):
    cell_pos_x = 0.5
    cell_pos_y = cell_pos_x

    num_w = len(sims[0])
    num_d = len(sims)
    data = np.zeros([num_w, num_d])

    for i,d  in enumerate(sorted(sims.keys())):
        for j, exp in enumerate(sims.get(d)):
            idx = exp.integrator.nPointsSpatial * cell_pos_x
            idy = exp.integrator.nPointsSpatial * cell_pos_y
            res = exp.singleCellTemporalResponse(cell, idx, idy)
            data[j][i] = np.mean(res)

    return data

R = area_summation(data, "relay")
I = area_summation(data, "interneuron")
G = area_summation(data, "ganglion")
C = area_summation(data, "cortical")

d = sorted(data.keys())
weights = np.linspace(2, 0.1, 10)
internpolation = "none"
cmap =cmaps.viridis

extent = [min(d), max(d), min(weights), max(weights)]
figure()
title("Ganglion")
imshow(G, extent = extent, origin="upper", interpolation=internpolation, aspect='auto', cmap =cmap)
xlabel(r"Spot diameter [deg]", fontsize= 16)
ylabel(r"$K_{ri}$",fontsize= 25)
colorbar()
if record : savefig(os.path.join(output_dir, "w_vs_d_ganglion.png"))

figure()
title("Interneuron")
imshow(I, extent = extent, origin="upper", interpolation=internpolation,aspect='auto', cmap =cmap)
xlabel(r"Spot diameter [deg]", fontsize= 16)
ylabel(r"$K_{ri}$",fontsize= 25)
colorbar()
if record : savefig(os.path.join(output_dir, "w_vs_d_interneuron.png"))

figure()
title("Relay")
imshow(R, extent = extent, origin="upper", interpolation=internpolation,aspect='auto', cmap =cmap)
xlabel(r"Spot diameter [deg]", fontsize= 16)
ylabel(r"$K_{ri}$",fontsize= 25)
colorbar()
if record : savefig(os.path.join(output_dir, "w_vs_d_relay.png"))


figure()
for i in range(1,len(weights),3):
    plot(d,G[i,:], label=r"$K_{ri}=$"+'{0:.2f}'.format(weights[i]))
legend()
title("Ganglion")
xlabel(r"Spot diameter [deg]", fontsize= 16)
ylabel("response",fontsize= 16)
if record : savefig(os.path.join(output_dir, "area_response_ganglion.png"))

figure()
for i in range(1,len(weights),3):
    plot(d,I[i,:], label=r"$K_{ri}=$"+'{0:.2f}'.format(weights[i]))
legend()
title("Interneuron")
xlabel(r"Spot diameter [deg]", fontsize= 16)
ylabel("response",fontsize= 16)
if record : savefig(os.path.join(output_dir, "area_response_interneuron.png"))


figure()
for i in range(1,len(weights),3):
    plot(d, R[i,:], label=r"$K_{ri}=$"+'{0:.2f}'.format(weights[i]))
legend()
title("Relay")
xlabel(r"Spot diameter [deg]", fontsize= 16)
ylabel("response",fontsize= 16)
if record : savefig(os.path.join(output_dir, "area_response_relay.png"))
