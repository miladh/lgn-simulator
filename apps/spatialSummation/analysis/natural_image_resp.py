#!/usr/bin/python
from pylab import*
import os, sys
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../../tools")))
from analysis.data_extractor import*
from analysis.pretty_plotting import*

#Analysis: ###########################################################################
def make_plot(cell_type, sims, cmap, save_fig=True):
    fig, axarr = plt.subplots(1, 3, figsize=(20,6))
    set_font()
    set_legend()
    for ax in axarr:
        spines_edge_color(ax)
        remove_ticks(ax)
        set_grid(ax, linewidth=0., linecolor='0.95')
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        spines_edge_color(ax, edges = {"top": "none", "bottom": "none",
                                       "right": "none", "left": "none"})



    w_rc = extract_unique_simulation_attrs(sims, "relay.Krc.w")
    w_ic = extract_unique_simulation_attrs(sims, "interneuron.Kic.w")
    print "w_rc=", w_rc
    print "w_ic=", w_ic

    sim_no_fb = simulation_extractor(sims, "relay.Krc.w", w_rc[0])[0]
    sim_with_fb = simulation_extractor(sims, "relay.Krc.w", w_rc[-1])[0]



    data1 =  np.copy(sim_no_fb.relay.resp()[0,:,:])
    data2 =  np.copy(sim_with_fb.relay.resp()[0,:,:])
    stim = sim_no_fb.stimulus.spatio_temporal()[0,:,:]

    vmax = max(abs(data2.min()), abs(data2.max()))
    vmin = -vmax

    axarr[0].imshow(stim,cmap="gray", aspect="auto", interpolation="none", origin="lower")
    axarr[0].set_title("Input",y=1.02)

    axarr[1].imshow(data1, vmin=vmin, vmax=vmax,
                   cmap=cmap, aspect="auto", interpolation="none", origin="lower")


    label = r"Without feedback"+"\n"+" $w_{\mathrm{RC}}=$"+'${0:.1f}$'.format(sim_no_fb.get_attribute("relay.Krc.w"))+r"$, w_{\mathrm{IC}}=$"+'${0:.1f}$'.format(sim_no_fb.get_attribute("interneuron.Kic.w"))
    axarr[1].set_title(label,y=0.99)


    label = r"With feedback"+"\n"+" $w_{\mathrm{RC}}=$"+'${0:.1f}$'.format(sim_with_fb.get_attribute("relay.Krc.w"))+r"$, w_{\mathrm{IC}}=$"+'${0:.1f}$'.format(sim_with_fb.get_attribute("interneuron.Kic.w"))
    im = axarr[2].imshow(data2, vmin=vmin, vmax=vmax,
                   cmap=cmap, aspect="auto", interpolation="none", origin="lower")


    plt.colorbar(im, ax=axarr.ravel().tolist(), orientation='vertical')
    axarr[2].set_title(label,y=0.99)

    ########################################################################################
    # plt.tight_layout()
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
    fig_name = "natural_image_fb"
    cmap="RdBu_r"

    sims = get_simulations(sims_path)
    Ns=sims[0].integrator.Ns
    Nt=sims[0].integrator.Nt
    s_points = sims[0].integrator.s_points[Ns/2:]
    k_points = sims[0].integrator.k_points[Ns/2:]
    rc = [Ns/2, Ns/2]
    #-----------------------------------------------------------------------------------

    for cell in cell_types:
        make_plot(cell, sims, cmap)
