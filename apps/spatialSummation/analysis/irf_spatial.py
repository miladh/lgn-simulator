#!/usr/bin/python
from pylab import*
import os, sys
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../../tools")))
from analysis.data_extractor import*

#Analysis: ###########################################################################
def extract_irfs(cell_type):
    global attr_a, attr_b
    irf_max = zeros([len(attr_a), len(attr_b)])
    irf_min = zeros([len(attr_a), len(attr_b)])
    irf_size = zeros([len(attr_a), len(attr_b)])

    irf_norm = norm_sim.get_attribute(cell_type).irf()[0, Ns/2,Ns/2:]
    for i, a in enumerate(attr_a):
        sims_a = simulation_extractor(sims, attr_a_name, a)
        for j, b in enumerate(attr_b):
            sim = simulation_extractor(sims_a, attr_b_name, b, return_as_list=False)
            irf = sim.get_attribute(cell_type).irf()[0, Ns/2, Ns/2:]

            irf_max[i,j]  = irf[0] / abs(irf_norm[0])
            irf_min[i,j]  = irf.min() / abs(irf_norm.min())
            irf_size[i,j] = s_points[argmin(irf)] / s_points[argmin(irf_norm)]


    return irf_max, irf_min, irf_size


def make_plot(irf_max, irf_min, irf_size, cell_type):
    import seaborn.apionly as sns
    sns.set_color_codes()
    from analysis.pretty_plotting import*

    fig, axarr = plt.subplots(1,3, figsize=(15,5), sharey='row')
    set_legend()
    set_font()
    for ax in axarr:
        spines_edge_color(ax)
        remove_ticks(ax)
        set_grid(ax, linewidth=0., linecolor='0.95')

    levels = levels = np.arange(-0.02, 1.0, 0.2)
    cmap="RdBu_r"
    extent = [abs(attr_b[0]), abs(attr_b[-1]),  abs(attr_a[0]), abs(attr_a[-1])]

    axarr[0].set_title("IRF center",y=1.02)
    im =axarr[0].imshow(irf_max, extent = extent,
                   cmap=cmap, aspect="auto",
                   interpolation="gaussian",
                   origin="lower")
    plt.colorbar(im, ax = axarr[0])

    axarr[1].set_title("IRF surround",y=1.02)
    im = axarr[1].imshow(irf_min, extent = extent,
               cmap=cmap, aspect="auto",
               interpolation="gaussian",
               origin="lower")
    plt.colorbar(im, ax = axarr[1])

    axarr[2].set_title("IRF size",y=1.02)
    im =axarr[2].imshow(irf_size, extent = extent,
               cmap=cmap, aspect="auto",
               interpolation="gaussian",
               origin="lower")
    plt.colorbar(im, ax = axarr[2])

    axarr[0].set_ylabel(ylabel)
    axarr[0].set_xlabel(xlabel)
    axarr[1].set_xlabel(xlabel)
    axarr[2].set_xlabel(xlabel)


    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, fig_name+cell_type+"_"+record_label+".pdf"))
    plt.show()


if __name__ == "__main__":
    import sumatra_tracking.io_manager as smt
    record_label = sys.argv[1:][-1]
    sims_path = sys.argv[1:][-2]
    output_dir = smt.get_output_dir(record_label)
    #-----------------------------------------------------------------------------------

    cell_type = ["relay", "interneuron"]
    attr_a_name = "interneuron.Kic.w"
    attr_b_name = "relay.Krc.w"

    xlabel ="$|w_{\mathrm{RC}}|$" #attr_b
    ylabel = "$|w_{\mathrm{IC}}|$" #attr_a
    fig_name= "irf_FB_"

    sims = get_simulations(sims_path)
    Ns=sims[0].integrator.Ns
    Nt=sims[0].integrator.Nt
    s_points = sims[0].integrator.s_points[Ns/2:]
    k_points = sims[0].integrator.k_points[Ns/2:]

    attr_a = extract_unique_simulation_attrs(sims, attr_a_name)
    attr_a = attr_a[argsort(abs(attr_a))]
    attr_a = attr_a[:]

    attr_b = extract_unique_simulation_attrs(sims, attr_b_name)
    attr_b = attr_b[argsort(abs(attr_b))]
    attr_b = attr_b[:]

    norm_sim = simulation_extractor(simulation_extractor(sims, attr_a_name, 0), attr_b_name, 0)[0]


    print "attr_a=", attr_a
    print "attr_b=", attr_b


    for cell in cell_type:
        make_plot(*extract_irfs(cell), cell_type=cell)
