#!/usr/bin/python
from pylab import*
import os, sys
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../../tools")))
from analysis.data_extractor import*
from analysis.pretty_plotting import*
import seaborn.apionly as sns
sns.set_color_codes()
#Analysis: ###########################################################################
def make_joint_plot(sim, save_fig=True):
    irf = sim.relay.irf()
    s_points = sim.integrator.s_points
    Ns=sim.integrator.Ns
    Nt=sim.integrator.Nt

    fig = plt.figure(figsize=(8,8))
    set_legend()
    set_font()

    gs = GridSpec(5,5)
    gs.update(wspace=0.025, hspace=0.05)

    ax_joint = fig.add_subplot(gs[2:5,0:3])
    ax_marg_x = fig.add_subplot(gs[0:2,0:3])
    ax_marg_y = fig.add_subplot(gs[2:5,3:5])
    axarr = [ax_marg_x,ax_marg_y]

    for ax in axarr:
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.spines['bottom'].set_position(('data',0))
        ax.yaxis.set_ticks_position('left')
        ax.spines['left'].set_position(('data',0))
        remove_ticks(ax)


    extent = [s_points.min(), s_points.max(), s_points.min(), s_points.max()]
    ax_joint.imshow(irf[0,:,:], extent = extent, cmap="RdBu_r", aspect="auto",
    interpolation="gaussian", origin="lower")
    ax_joint.set_xlim([-2, 2])
    ax_joint.set_ylim([-2, 2])
    spines_edge_color(ax_joint)
    remove_ticks(ax_joint)
    set_grid(ax_joint, linewidth=0., linecolor='0.95')
    spines_edge_color(ax_joint, edges = {"top": "none", "bottom": "none",
                                   "right": "none", "left": "none"})


    ax_marg_x.plot(s_points, irf[0,Ns/2,:])
    ax_marg_x.set_xlim([-2, 2])
    # ax_marg_x.set_ylim([-0.1, 0.7])

    ax_marg_y.plot( irf[0,:,Ns/2],s_points)
    # ax_marg_y.set_xlim([-3, 3])
    ax_marg_y.set_ylim([-2, 2])


    # Turn off tick labels on marginals
    plt.setp(ax_marg_x.get_xticklabels(), visible=False)
    plt.setp(ax_marg_x.get_yticklabels(), visible=False)
    plt.setp(ax_marg_y.get_xticklabels(), visible=False)
    plt.setp(ax_marg_y.get_yticklabels(), visible=False)

    ax_joint.set_xlabel('$x (^\circ)$')
    ax_joint.set_ylabel('$y (^\circ)$')
    # ax_marg_y.set_xlabel('$W(x=0, y)$')
    # ax_marg_x.set_ylabel('$W(x, y=0)$')


    irf_max = max(irf[0,Ns/2,:])
    irf_min = min(irf[0,Ns/2,:])
    irf_size = argmin(irf[0,Ns/2,:])


    ax_marg_x.scatter(0,irf_max, s=80, marker="s",facecolor="r" ,
    label="Center excitation", zorder=3)
    ax_marg_x.scatter(s_points[-irf_size],irf_min, facecolor="r", s=80,
    label="Surround inhibiton", zorder=3)

    ax_marg_x.scatter(s_points[-irf_size], 0.03, s=80, marker="v",facecolor="r" ,
    label="Size", zorder=3)

    ax_marg_x.plot([s_points[-irf_size],s_points[-irf_size]],[0,irf_min],
    color ='r', linewidth=1.5, linestyle="--")
    ax_marg_x.plot([0, s_points[-irf_size]],[irf_min,irf_min], color ='r',
     linewidth=1.5, linestyle="--")

    ax_marg_x.legend(loc=(0.7,0.7), fontsize=14)

    if save_fig: fig.savefig(os.path.join(output_dir, fig_name+cell_type+"_"+record_label+".pdf"))
    if save_fig: fig.savefig(os.path.join(output_dir, fig_name+cell_type+"_"+record_label+".png"))
    plt.show()


if __name__ == "__main__":
    import sumatra_tracking.io_manager as smt
    record_label = sys.argv[1:][-1]
    sims_path = sys.argv[1:][-2]
    output_dir = smt.get_output_dir(record_label)

    #-----------------------------------------------------------------------------------
    cell_type = ["relay"]
    fig_name= "edog_spatial_irf_"
    sims = get_simulations(sims_path)
    Ns=sims[0].integrator.Ns
    Nt=sims[0].integrator.Nt
    s_points = sims[0].integrator.s_points[Ns/2:]
    sim = sims[0]
    #-----------------------------------------------------------------------------------
    make_joint_plot(sim)
