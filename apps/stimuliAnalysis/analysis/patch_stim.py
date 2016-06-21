#!/usr/bin/python
import os, sys
from argparse import ArgumentParser
from pylab import*


current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../../tools")))

import sumatra_tracking.io_manager as smt
import analysis.colormaps as cmaps
from analysis.tuning_analysis import*
from analysis.data_extractor import*
from analysis.pretty_plotting import*
from analysis.data_extractor import*

parser = ArgumentParser()
parser.add_argument("sim_ids", help = "simulation ids")
parser.add_argument("record", help = "record results", type = int)
parser.add_argument("run_id", help = "sumatra_label")
args = parser.parse_args()
sim_ids = args.sim_ids
record = args.record
run_id = args.run_id


output_dir = None
if(record):
    output_dir = smt.get_output_dir(sim_ids, run_id)

sims = get_simulations(sim_ids)
# Analysis: --------------------------------------------------------------------
Ns=sims[-1].integrator.Ns
Nt=sims[-1].integrator.Nt
for sim in sims:
    sim.stimulus.read_property()

diameters = extract_unique_simulation_attrs(sims, "stimulus.mask_size")
diameters = diameters[argsort(diameters)]
k_points = sims[0].integrator.k_points
print diameters

fig = plt.figure()

ax = fig.add_subplot(111)
spines_edge_color(ax)
remove_ticks(ax)
set_grid(ax)
set_font()
set_legend()


n=[0,4,6]
for i, d in enumerate(diameters[2:-1:10]):
    label=r"$d=$"+'{0:.2f}'.format(d)

    sim = simulation_extractor(sims, "stimulus.mask_size", d)[0]
    ax.plot(sim.integrator.k_points,
    sim.stimulus.fourier_transform[0, Ns/2,:]/max(sim.stimulus.fourier_transform[0, Ns/2, :]),
    label = label, color=colormap(n[i]))

ax.set_title("$\widetilde{S}(k_x, k_y=0, w=0; d)$")
ax.set_xlabel("$k_x$")
ax.set_ylabel("$\widetilde{S}$")
ax.set_xlim([-30, 30])
legend()
if record : fig.savefig(os.path.join(output_dir, "stim_ft.png"))

# 2d: --------------------------------------------------------------------
S_ft = np.zeros([len(diameters), Ns])

for i, d in enumerate(diameters):
    sim = simulation_extractor(sims, "stimulus.mask_size", d)[0]
    S_ft[i,:] = sim.stimulus.fourier_transform[0, Ns/2]/max(sim.stimulus.fourier_transform[0, Ns/2])

fig = plt.figure()
extent =[k_points.min(), k_points.max(), diameters.min(), diameters.max()]
internpolation = "gaussian"
imshow(S_ft, extent=extent, origin="lower", aspect='auto', interpolation=internpolation, cmap =cmaps.viridis )
title("$\widetilde{S}(k_x, k_y=0, w=0; d)$")
ylabel("Spot diameter [deg]", fontsize= 16)
xlabel("$k_x$",fontsize= 25)
colorbar()
if record : fig.savefig(os.path.join(output_dir, "stim_ft_vs_d.png"))
