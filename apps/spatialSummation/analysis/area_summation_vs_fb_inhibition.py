#!/usr/bin/python
import os, sys
import sumatra_tracking.io_manager as smt

current_path = os.path.dirname(os.path.realpath(__file__))
lib_path = [os.path.abspath(os.path.join(current_path,"../../../tools/sumatraTracking")),
           os.path.abspath(os.path.join(current_path,"../../../tools/analysis"))]
[sys.path.append(path) for path in lib_path]

options = sys.argv[1:]
record_label = options[-1]
sims_path = options[-2]
output_dir = smt.get_output_dir(record_label)



# Analysis: --------------------------------------------------------------------
import seaborn.apionly as sns
sns.set_color_codes()
from analysis.data_extractor import*
from analysis.pretty_plotting import*

sims = get_simulations(sims_path)
s_points = sims[0].integrator.s_points
Ns=sims[0].integrator.Ns
Nt=sims[0].integrator.Nt

a_rg = extract_unique_simulation_attrs(sims, "relay.Krg.spatial.a", return_as_array=False)

weights = extract_unique_simulation_attrs(sims, "relay.Krc.w")
weights = weights[argsort(weights)]
weights = weights[:-4]

widths = extract_unique_simulation_attrs(sims, "relay.Krc.spatial.a")
widths = widths[argsort(widths)]
widths = widths[:]


print "w_rc=", weights
print "a_rc=", widths
print "a_rg=", a_rg


# Tests:
print "Kri.w=", extract_unique_simulation_attrs(sims, "relay.Kri.w")
print "Kri.a=", extract_unique_simulation_attrs(sims, "relay.Kri.spatial.a")
print "Kic.w=" ,extract_unique_simulation_attrs(sims, "interneuron.Kic.w")
print "Kic.a=" ,extract_unique_simulation_attrs(sims, "interneuron.Kic.spatial.a")
print "Kig.w=",extract_unique_simulation_attrs(sims, "interneuron.Kig.w")
print "Kig.a=",extract_unique_simulation_attrs(sims, "interneuron.Kig.spatial.a")

def make_plot(cell_type):
    rc_max = zeros([len(widths), len(weights)])
    rc_min = zeros([len(widths), len(weights)])
    rc_size = zeros([len(widths), len(weights)])

    no_inhib = simulation_extractor(sims, "relay.Krc.w", 0)[0].get_attribute(cell_type).irf()[0, Ns/2,:]
    for i, w in enumerate(weights):
        sim_w = simulation_extractor(sims, "relay.Krc.w", w)
        #########################################################################################################################
        for j, a in enumerate(widths[:]):
            sim = simulation_extractor(sim_w, "relay.Krc.spatial.a", a, return_as_list=False)
            data = sim.get_attribute(cell_type).irf()[0, Ns/2,:]
            rc_max[j, i] = data[Ns/2]/no_inhib.max()

            if(rc_max[j, i]<0):
                rc_min[j, i] = data[Ns/2+1:].max()/abs(no_inhib.min())
                rc_size[j, i] = abs(s_points[where(data==data[Ns/2+1:].max())][0])/abs(s_points[where(no_inhib==no_inhib.min())][0])

            else:
                rc_min[j, i] = data.min()/abs(no_inhib.min())
                rc_size[j, i] = abs(s_points[where(data==data.min())][0])/abs(s_points[where(no_inhib==no_inhib.min())][0])


    fig, axarr = plt.subplots(1,3, figsize=(15,5), sharey='row')
    set_legend()
    set_font()
    for ax in axarr:
        spines_edge_color(ax)
        remove_ticks(ax)
        set_grid(ax, linewidth=0., linecolor='0.95')

    levels = levels = np.arange(-0.02, 1.0, 0.2)
    cmap="RdBu_r"
    extent = [ 0, 1, widths.min(), widths.max()]

    axarr[0].set_title("Response center",y=1.02)
    im =axarr[0].imshow(rc_max, extent = extent,
                   #vmin=vmin, vmax=vmax,
                   cmap=cmap, aspect="auto",
                   interpolation="gaussian",
                   origin="lower")
    plt.colorbar(im, ax = axarr[0])

    axarr[1].set_title("Response surround",y=1.02)
    im = axarr[1].imshow(rc_min, extent = extent,
               #vmin=vmin, vmax=vmax,
               cmap=cmap, aspect="auto",
               interpolation="gaussian",
               origin="lower")
    plt.colorbar(im, ax = axarr[1])

    axarr[2].set_title("RF size",y=1.02)
    im =axarr[2].imshow(rc_size, extent = extent,
               #vmin=vmin, vmax=vmax,
               cmap=cmap, aspect="auto",
               interpolation="gaussian",
               origin="lower")
    plt.colorbar(im, ax = axarr[2])

    axarr[0].set_ylabel("$a_{\mathrm{RC}}$")
    axarr[0].set_xlabel("$|w_{\mathrm{RC}}|$")
    axarr[1].set_xlabel("$|w_{\mathrm{RC}}|$")
    axarr[2].set_xlabel("$|w_{\mathrm{RC}}|$")


    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "irf_noFB_"+cell_type+"_"+run_id+".pdf"))

    plt.show()



cell_types = ["relay", "interneuron", "ganglion"]

for cell in cell_types[:1]:
    make_plot(cell)
