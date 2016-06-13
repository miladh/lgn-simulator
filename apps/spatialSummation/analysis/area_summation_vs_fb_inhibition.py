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

sims = get_simulations(sim_ids)

output_dir = None
if(record):
    output_dir = smt.get_output_dir(sim_ids, run_id)

# Analysis: --------------------------------------------------------------------
cell_types= {"relay": "relay.Krc.w",
            "interneuron": "interneuron.Kic.w",
            "cortical": "relay.Krc.w", }

for cell, attr in cell_types.iteritems():
    print cell, attr
    weights = extract_unique_simulation_attrs(sims, attr)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    spines_edge_color(ax)
    remove_ticks(ax)
    set_grid(ax)
    set_font()
    set_legend()

    n=0
    for w in weights:
        n+=1
        sims_ext = simulation_extractor(sims, attr, w)
        data=average_response_vs_attr(sims_ext, "stimulus.mask_size")
        sorted_indices = argsort(data["stimulus.mask_size"])
        mask_size = array(data["stimulus.mask_size"])[sorted_indices]
        resp = array(data[cell])[sorted_indices]
        ax.plot(mask_size, resp, "-^", color=colormap(n),
                label=r"$K_{rc}=$"+'{0:.2f}'.format(w))

    ax.set_title(cell)
    ax.set_xlabel("diameter[deg]")
    ax.set_ylabel("response[spikes/sec]")
    if record : fig.savefig(os.path.join(output_dir, "area_response"+"_"+cell+".png"))
    legend()
