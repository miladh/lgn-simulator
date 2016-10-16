#!/usr/bin/python
from pylab import*
import os, sys
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../../tools")))
from analysis.data_extractor import*
from analysis.pretty_plotting import*
import seaborn.apionly as sns
sns.set_color_codes()
import matplotlib.colors as colors

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

#Analysis: ###########################################################################
def xy_plot(sim, cell, x_lim,  times, num_levels=[5,5], norm_id=0, save_fig=True):
    Ns=sim.integrator.Ns
    Nt=sim.integrator.Nt
    s_points = sim.integrator.s_points
    t_points = sim.integrator.t_points
    extent = [s_points.min(), s_points.max(), s_points.min(), s_points.max()]

    fig, axarr = plt.subplots(2, len(times), figsize=(3*len(times), 6), sharey="row")

    set_legend()
    set_font()
    cmap="RdBu_r"
    for ax in axarr[0,:]:
        spines_edge_color(ax)
        remove_ticks(ax)
        set_grid(ax, linewidth=0., linecolor='0.95')
        spines_edge_color(ax, edges = {"top": "none", "bottom": "none",
        "right": "none", "left": "none"})

    for ax in axarr[1,:]:
        remove_ticks(ax)
        move_spines(ax)
        remove_ticklabels(ax)
        spines_edge_color(ax, edges = {"top": "none", "bottom": "none",
                                          "right": "none", "left": "none"})

    vmax = (sim.get_attribute(cell).irf()[times[norm_id],:,:]).max()
    vmin = (sim.get_attribute(cell).irf()[times[norm_id],:,:]).min()
    levels = unique(append(linspace(vmin,0, num_levels[0]), linspace(0., vmax, num_levels[1])))
    print levels
    for i, t in enumerate(times):
        ax = axarr[0,i]
        label = r"$t=$"+'${0:.0f}$'.format(t_points[times[i]])+"$\mathrm{ms}$"
        irf = sim.get_attribute(cell).irf()[t,:,:]
        X, Y = np.mgrid[ s_points.min():s_points.max():complex(0, Ns), s_points.min():s_points.max():complex(0, Ns)]
        im =ax.pcolormesh(Y, X, irf, norm=MidpointNormalize(midpoint=0.),
            cmap=cmap, vmin=vmin, vmax=vmax)
        ax.contour( irf, extent= extent,levels=levels,
                    colors='k',linewidths=0.4,
                    aspect="auto")
        ax.set_xlim([-x_lim, x_lim])
        ax.set_ylim([-x_lim, x_lim])
        ax.set_xlabel("$x(^\circ)$", fontsize=18)
        ax.set_title(label,y=1, fontsize=18)

        ax = axarr[1,i]
        ax.plot(s_points, irf[Ns/2,:], "k", linewidth=1.5)
        ax.set_xlim([-x_lim, x_lim])
        ax.fill_between(s_points, 0, irf[Ns/2,:], where=irf[Ns/2,:] < 0, facecolor='#5B9FCF', alpha=0.7)
        ax.fill_between(s_points, 0, irf[Ns/2,:],where=irf[Ns/2,:] >= 0,  facecolor='#E65D6D', alpha=0.7)



    axarr[0,0].set_ylabel("$y(^\circ)$", fontsize=18)
    plt.tight_layout()
    if save_fig: fig.savefig(os.path.join(output_dir, fig_name+cell_type+"_"+record_label+".pdf"))
    if save_fig: fig.savefig(os.path.join(output_dir, fig_name+cell_type+"_"+record_label+".png"))
    plt.show()

#------------------------------------------------------------------------------------------------------------
def xt_plot(sims, cell_type, x_lim, t_lim,  levels, save_fig=True):
    Ns=sims[0].integrator.Ns
    Nt=sims[0].integrator.Nt
    s_points = sims[0].integrator.s_points
    t_points = sims[0].integrator.t_points
    extent = [s_points.min(), s_points.max(), t_points.min(), t_points.max()]
    weights = extract_unique_simulation_attrs(sims, "relay.Krc.w")
    weights = weights[argsort(weights)]

    fig, axarr = plt.subplots(1, len(sims), figsize=(4*len(sims), 5), sharey="row")

    set_legend()
    set_font()
    cmap="RdBu_r"

    if len(sims)==1: axarr=[axarr]
    for ax in axarr:
        spines_edge_color(ax)
        remove_ticks(ax)
        set_grid(ax, linewidth=0., linecolor='0.95')
        spines_edge_color(ax, edges = {"top": "none", "bottom": "none",
        "right": "none", "left": "none"})


    ################################################################################
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    for ax in axarr:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

    for i in range(len(fig.axes), len(axarr),-1):
        fig.delaxes(fig.axes[i-1])

    divider = make_axes_locatable(axarr[-1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    ###############################################################################


    for i, w in enumerate(weights):
        sim = simulation_extractor(sims, "relay.Krc.w", w, return_as_list=False)
        ax = axarr[i]
        irf = sim.get_attribute(cell_type).irf()[:,Ns/2,:]
        X, Y = np.mgrid[ t_points.min():t_points.max():complex(0, Nt), s_points.min():s_points.max():complex(0, Ns)]
        im =ax.pcolormesh(Y, X, irf, norm=MidpointNormalize(midpoint=0.), cmap=cmap,                        vmin=-0.1, vmax=0.25)
        ax.contour( irf, extent= extent,
                    colors='k',linewidths=0.4,
                    aspect="auto",
                    levels = levels)

        ax.set_xlim([-x_lim, x_lim])
        ax.set_ylim([0, t_lim])
        ax.set_xlabel("$x(^\circ)$", fontsize=18)

        label = r"$w_{\mathrm{RCR}}=$"+'${0:.2f}$'.format(sim.get_attribute("relay.Krc.w"))
        ax.set_title(label,y=1, fontsize=20)



    plt.colorbar(im, cax = cax, ax = axarr[-1])
    axarr[0].set_ylabel(r"$\tau\;(\mathrm{ms})$", fontsize=18)
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
