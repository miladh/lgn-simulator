#!/usr/bin/python
from pylab import*
import os, sys
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../../tools")))
from analysis.data_extractor import*
from analysis.tuning_analysis import*
from analysis.pretty_plotting import*
import seaborn.apionly as sns
sns.set_color_codes()


#Analysis: ###########################################################################
def make_plot(cell_type, resp, attr_a, attr_b, freqs, diameter, save_fig=True):
    fig, ax = plt.subplots(1, 1, figsize=(6,6),  sharex='col')
    set_font()
    set_legend()
    spines_edge_color(ax)
    remove_ticks(ax)
    set_grid(ax)

    for wi, wr in zip(attr_a[[0, 1, 3]], attr_b[[0, 1, 3]]):
        i = where(attr_a==wi)[0][0]
        label = r"$w_{\mathrm{RCR}}=$"+'${0:.2f}$'.format(wr)+r"$, w_{\mathrm{ICR}}=$"+'${0:.1f}$'.format(wi)
        ax.plot(freqs, resp[cell_type][i,:], "-", label=label)

    ax.set_ylabel("Response", fontsize=20)
    ax.set_xlabel("Patch-grating wave vector $k_\mathrm{pg} (1/^\circ)$", fontsize=20)
    ax.set_title("$\mathrm{Patch\;size}=$"+'${0:.2f}^\circ$'.format(diameter), fontsize=20)
    ax.legend()
    ax.set_xlim([0, 10])


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
    cell_types = ["relay"]
    attr_a_name = "interneuron.Kic.w"
    attr_b_name = "relay.Krc.w"
    fig_name= "freq_tuning_fb_weights_small_d"

    sims = get_simulations(sims_path)
    Ns=sims[0].integrator.Ns
    Nt=sims[0].integrator.Nt
    s_points = sims[0].integrator.s_points[Ns/2:]
    k_points = sims[0].integrator.k_points[Ns/2:]
    rc = [Ns/2, Ns/2]

    attr_a = extract_unique_simulation_attrs(sims, attr_a_name)
    attr_b = extract_unique_simulation_attrs(sims, attr_b_name)
    attr_a2, freqs, resp = resp_vs_attrA_vs_attrB(sims, attr_a_name, "stimulus.spatial_freq", rc=rc)
    if(dot(attr_a - attr_a2,attr_a - attr_a2)!=0 ):
        raise ValueError('attra and attra2 are different:' + attr_a, attr_a2)
    else:
       print "diff=", attr_a - attr_a2

    d = extract_unique_simulation_attrs(sims, "stimulus.mask_size")
    if len(d)>1:raise IndexError(d)
    d=d[0]
    print d

    #-----------------------------------------------------------------------------------

    for cell in cell_types:
        make_plot(cell, resp, attr_a, attr_b, freqs, d)
