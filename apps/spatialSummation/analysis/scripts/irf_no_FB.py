
# coding: utf-8

# ###Analyzing the impulse response function:

# In[1]:


from pylab import*
import seaborn as sns
import os, sys

current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../../../tools")))

import sumatra_tracking.io_manager as smt
from analysis.data_extractor import*
from analysis.pretty_plotting import*

run_id = "20160815-180548"
sim_ids = "/media/milad/scratch/lgn-simulator/simulations/spatialSummation/" + run_id
sims=get_simulations(sim_ids)


# In[2]:

s_points = sims[0].integrator.s_points
Ns=sims[0].integrator.Ns
Nt=sims[0].integrator.Nt

a_rg = extract_unique_simulation_attrs(sims, "relay.Krg.spatial.a", return_as_array=False)

weights = extract_unique_simulation_attrs(sims, "relay.Kri.w")
weights = weights[argsort(-weights)]
weights = weights[1:]

widths = extract_unique_simulation_attrs(sims, "relay.Kri.spatial.a")
widths = widths[argsort(widths)]
widths = widths[:]


print "w_ri=", weights
print "a_ri=", widths
print "a_rg=", a_rg


# In[3]:

def make_plot(cell_type):
    fig, axarr = plt.subplots(2,len(weights), figsize=(12,8),  sharex='col', sharey="row")
    set_legend()
    set_font()
    for ax_row in axarr:
        for ax in ax_row:
            spines_edge_color(ax)
            remove_ticks(ax)
            set_grid(ax)


    no_inhib = simulation_extractor(sims, "relay.Kri.w", 0)[0].get_attribute(cell_type).irf()[0, Ns/2,:]
    for i, w in enumerate(weights):
        sim_w = simulation_extractor(sims, "relay.Kri.w", w)
        #########################################################################################################################
        ax = axarr[0, i]
        ax.set_title("$w_\mathrm{RI}="+"{0:.2f}$".format(w), y=1.02, fontsize=18)
        ax.plot(s_points, no_inhib,  "--k", label ="No inhibition")
        for a in widths:
            label="$a_\mathrm{RI}/a_\mathrm{RG}=$"+"${0:.1f}$".format(a/a_rg)
            sim = simulation_extractor(sim_w, "relay.Kri.spatial.a", a, return_as_list=False)
            ax.plot(sim.integrator.s_points, sim.get_attribute(cell_type).irf()[0, Ns/2,:], label = label)
       ###########################################################################################################################
        ax = axarr[1, i]
        ax.set_xlim(0,2)
        ax.set_xlabel("$x\;(^\circ)$",  fontsize=14)
        ax.plot(s_points, no_inhib/no_inhib.max(),  "--k", label ="No inhibition")
        for a in widths:
            sim = simulation_extractor(sim_w, "relay.Kri.spatial.a", a, return_as_list=False)
            ax.plot(sim.integrator.s_points,
            sim.get_attribute(cell_type).irf()[0, Ns/2,:]/max(sim.get_attribute(cell_type).irf()[0, Ns/2,:]),
            label = label)


    axarr[0,0].legend(fontsize=14)
    axarr[0,0].set_ylabel("$W_\mathrm{R}(x, y=0)$",  fontsize=14)
    axarr[1,0].set_ylabel("$W_\mathrm{R}(x, y=0)$ (Normalized)",  fontsize=14)
    plt.tight_layout()
    #plt.savefig("/home/milad/Dropbox/projects/lgn/paper/images/results/irf_noFB_"+cell_type+"_"+run_id+".pdf")
    plt.show()




# In[5]:

cell_types = ["relay", "interneuron", "ganglion"]

for cell in cell_types[:1]:
    make_plot(cell)


# In[4]:
