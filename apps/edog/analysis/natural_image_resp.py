#!/usr/bin/python
from pylab import*
import os, sys
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../../tools")))
from analysis.data_extractor import*
from analysis.pretty_plotting import*

import matplotlib.colors as colors
class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

#Analysis: ###########################################################################
def make_plot(cell_type, sims, cmap, save_fig=True):
    fig, axarr = plt.subplots(1, 3, figsize=(14,5))
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
    w_rc_c = extract_unique_simulation_attrs(sims, "relay.Krc.spatial.c")
    print "w_rc=", w_rc
    print "w_ic=", w_rc_c

    sim_no_fb = simulation_extractor(sims, "relay.Krc.w", w_rc[0])[0]
    sim_with_fb = simulation_extractor(sims, "relay.Krc.w", w_rc[-1])[0]



    data1 =  np.copy(sim_no_fb.relay.resp()[0,:-50,:])
    data2 =  np.copy(sim_with_fb.relay.resp()[0,:-50,:])
    stim = sim_no_fb.stimulus.spatio_temporal()[0,:-50,:]


    vmax = max(abs(data2.min()), abs(data2.max()))
    vmin = -vmax

    im1=axarr[0].imshow(stim,cmap="gray", aspect="auto", interpolation="none", origin="lower")
    axarr[0].set_title("Input",y=1.02, fontsize=20)

    im2=axarr[1].imshow(data1, vmin=vmin, vmax=vmax,
                   cmap=cmap, aspect="auto", interpolation="none", origin="lower")


    label = r"Without feedback"+"\n"+" $w_{\mathrm{RCR}}=$"+'${0:.1f}$'.format(sim_no_fb.get_attribute("relay.Krc.w"))
    axarr[1].set_title(label,y=1, fontsize=20)


    label = r"With feedback"+"\n"+" $w_{\mathrm{RCR}}=$"+'${0:.1f}$'.format(sim_with_fb.get_attribute("relay.Krc.w"))
    im3 = axarr[2].imshow(data2, vmin=vmin, vmax=vmax,
                   cmap=cmap, aspect="auto", interpolation="none", origin="lower")

    axarr[2].set_title(label,y=1, fontsize=20)
    ################################################################################
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider1 = make_axes_locatable(axarr[0])
    cax1 = divider1.append_axes("right", size="5%", pad=0.05)

    divider2 = make_axes_locatable(axarr[1])
    cax2 = divider2.append_axes("right", size="5%", pad=0.05)

    divider3 = make_axes_locatable(axarr[2])
    cax3 = divider3.append_axes("right", size="5%", pad=0.05)

    fig.delaxes(fig.axes[4])
    fig.delaxes(fig.axes[3])
    ###############################################################################

    cbar3 = fig.colorbar(im3, cax = cax3, ticks=[vmin,0,  vmax])
    cbar3.set_ticklabels(['Low', '0', 'High'])
    cbar3.ax.tick_params(labelsize=18)

    ########################################################################################
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
    attr_a_name = "relay.Kic.w"
    attr_b_name = "relay.Krc.w"
    fig_name = "natural_image_fb"
    cmap="RdBu_r"
    #cmap="Blues"
    #cmap="gray_r"


    sims = get_simulations(sims_path)
    Ns=sims[0].integrator.Ns
    Nt=sims[0].integrator.Nt
    s_points = sims[0].integrator.s_points[Ns/2:]
    k_points = sims[0].integrator.k_points[Ns/2:]
    rc = [Ns/2, Ns/2]
    #-----------------------------------------------------------------------------------

    for cell in cell_types:
        make_plot(cell, sims, cmap)
