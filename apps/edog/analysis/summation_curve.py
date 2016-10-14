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
def make_plot(cell_type, resp, attr_a, attr_b, diameter, save_fig=True):
    fig, axarr = plt.subplots(1, 3, figsize=(15,6))
    set_font()
    set_legend()
    for ax in axarr:
        spines_edge_color(ax)
        remove_ticks(ax)
        set_grid(ax)


    ax = axarr[0]
    for wr in attr_a[:-1]:
        i = where(attr_a==wr)[0][0]
        label = r"$w_{\mathrm{RCR}}=$"+'${0:.2f}$'.format(wr)
        ax.plot(diameter, resp[cell_type][i,:], "-", label=label)

    ax.set_xlabel("Diameter($^\circ$)")
    ax.set_ylabel("Response")
    ax.set_title("Area-response curve",y=1.02)
    ax.set_xlim([0., 10])
    # ax.set_ylim([0, 0.3])
    ax.legend()

    #########################################################################################
    ax = axarr[1]
    spines_edge_color(ax)
    remove_ticks(ax)
    ax.yaxis.set_tick_params(size=0)
    ax.xaxis.set_tick_params(size=0)

    d_max = zeros(len(attr_a))
    d_min = zeros(len(attr_a))
    for i, w in enumerate(attr_a[:]):
        index_max = argmax(resp[cell_type][i,0:])
        index_min = argmin(resp[cell_type][i,index_max:])
        d_max[i] = diameter[index_max]
        d_min[i] = diameter[index_min]

    ax.plot(attr_b, d_max, "-", label="Center")
    ax.plot(attr_b, d_min, "-", label="Surround")
    ax.set_ylabel("Diameter($^\circ$)")
    ax.set_xlabel("$w_{\mathrm{RCR}}$", fontsize=20)
    ax.tick_params(direction='out', pad=7)
    ax.legend()

    #########################################################################################
    ax = axarr[2]

    spines_edge_color(ax)
    remove_ticks(ax)
    ax.yaxis.set_tick_params(size=0)
    ax.xaxis.set_tick_params(size=0)

    supp = zeros(len(attr_a))
    for i, w in enumerate(attr_a):
        supp[i] = 1 - resp[cell_type][i,-1]/resp[cell_type][i,:].max()

    ax.plot(attr_b[:], supp[:], "-", label="Suppression")
    ax.set_xlabel("$w_{\mathrm{RCR}}$", fontsize=20)
    ax.set_ylabel("Suppresion index(%)")
    ax.set_ylim([0.75, 0.95])
    ax.tick_params(direction='out', pad=7)
    ax.set_yticks(arange(0.75, 0.95, 0.025))
    ax.set_yticks(arange(0.75, 0.975, 0.025), minor=True)
    # ax.grid(which='minor',color="w", linestyle="-", linewidth=1.3, zorder = 0)
    #########################################################################################
    plt.tight_layout()
    fig.savefig("area_summation_fb_weights.pdf")
    # if save_fig: fig.savefig(os.path.join(output_dir, fig_name+cell_type+"_"+record_label+".pdf"))
    # if save_fig: fig.savefig(os.path.join(output_dir, fig_name+cell_type+"_"+record_label+".png"))
    plt.show()



if __name__ == "__main__":
    import sumatra_tracking.io_manager as smt
    record_label = sys.argv[1:][-1]
    sims_path = sys.argv[1:][-2]
    output_dir = smt.get_output_dir(record_label)

    #-----------------------------------------------------------------------------------
    cell_types = ["relay"]
    attr_a_name = "relay.Krc.w"
    attr_b_name = "relay.Krc.w"
    diameters =  "stimulus.mask_size"
    fig_name= "area_summation_fb_weights_"

    sims = get_simulations(sims_path)
    Ns=sims[0].integrator.Ns
    Nt=sims[0].integrator.Nt
    s_points = sims[0].integrator.s_points[Ns/2:]
    k_points = sims[0].integrator.k_points[Ns/2:]
    rc = [Ns/2, Ns/2]

    attr_a = extract_unique_simulation_attrs(sims, attr_a_name)
    attr_b = extract_unique_simulation_attrs(sims, attr_b_name)
    attr_a2, diameters, resp = resp_vs_attrA_vs_attrB(sims, attr_a_name, diameters, rc=rc)
    if(dot(attr_a - attr_a2,attr_a - attr_a2)!=0 ):
        raise ValueError('attra and attr2 are different:' + attr_a, attr_a2)
    else:
       print "diff=", attr_a - attr_a2

    #-----------------------------------------------------------------------------------

    for cell in cell_types:
        make_plot(cell, resp, attr_a, attr_b, diameters)
