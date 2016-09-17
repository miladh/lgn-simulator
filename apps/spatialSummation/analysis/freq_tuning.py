#!/usr/bin/python
from pylab import*
import os, sys
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../../tools")))
from analysis.data_extractor import*
from analysis.tuning_analysis import*
from analysis.pretty_plotting import*

#Analysis: ###########################################################################
def make_plot(cell_type, resp, attr_a, attr_b, freqs, diameter, save_fig=True):
    fig, axarr = plt.subplots(2, 1, figsize=(5,10),  sharex='col')
    set_font()
    set_legend()
    for ax in axarr:
        spines_edge_color(ax)
        remove_ticks(ax)
        set_grid(ax)

    for wi, wr in zip(attr_a[[0, 1, 3]], attr_b[[0, 1, 3]]):
        i = where(attr_a==wi)[0][0]
        label = r"$w_{\mathrm{RC}}=$"+'${0:.2f}$'.format(wr)+r"$, w_{\mathrm{IC}}=$"+'${0:.1f}$'.format(wi)
        axarr[0].plot(freqs, resp[cell_type][i,:], "-", label=label)
        axarr[1].plot(freqs, resp[cell_type][i,:]/resp[cell_type][i,:].max() , "-", label=label)

    axarr[0].set_ylabel("Response(spikes/s)")
    axarr[0].set_title("$d=$"+'${0:.2f}$'.format(diameter))
    axarr[0].legend()
    axarr[1].set_ylabel("Normalized response(%)")
    #axarr[1].set_xscale('log')
    #axarr[1].set_xticks([0.5, 1, 2, 4, 8])
    #axarr[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axarr[1].set_xlim([0, 10])
    #axarr[1].legend()


    #########################################################################################
    plt.tight_layout()
    if save_fig: fig.savefig(os.path.join(output_dir, fig_name+cell_type+"_"+record_label+".pdf"))
    if save_fig: fig.savefig(os.path.join(output_dir, fig_name+cell_type+"_"+record_label+".png"))
    plt.show()



if __name__ == "__main__":
    import sumatra_tracking.io_manager as smt
    record_label = sys.argv[1:][-1]
    sims_path = sys.argv[1:][-2]
    output_dir = smt.get_output_dir(record_label)

    #-----------------------------------------------------------------------------------


    #-----------------------------------------------------------------------------------

    for cell in cell_types:
        make_plot(cell, resp, attr_a, attr_b, freqs, d)
