#!/usr/bin/python
from pylab import*
import os, sys
import h5py
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../tools")))
from analysis.data_extractor import*
from analysis.plotting_tools import*



data_file_path ="/home/milad/Dropbox/projects/lgn/code/build-lgnSimulator-Desktop_Qt_5_7_0_GCC_64bit-Debug/apps/default/out.h5"
config_file_path ="defaultConfig.yaml"
data_file = h5py.File(data_file_path, "r")

sim = Simulation(config_file_path, data_file)

Ns = sim.integrator.Ns
Nt = sim.integrator.Nt
print "Ns:", Ns, "Nt:", Nt


print "ds:", sim.integrator.ds
print "ds:", sim.integrator.dk
print "dt:", sim.integrator.dt
print "dt:", sim.integrator.dw

stim = sim.stimulus.spatio_temporal()
print stim.shape
plot(sim.integrator.t_points, stim[:, Ns/2,Ns/2])
show()

relay_resp = sim.relay.resp()
print relay_resp.shape
plot(sim.integrator.t_points, relay_resp[:, Ns/2,Ns/2])
show()

#Animation-------------------------------------------------------
S = {"type" : "Stimulus",
        "value" : sim.stimulus.spatio_temporal(),
        "t_points" : sim.integrator.t_points,
        "spatial_vec" : sim.integrator.s_points
        }

Rg = {"type" : "Ganglion",
            "value" : sim.ganglion.resp(),
            "t_points" : sim.integrator.t_points,
            "spatial_vec" : sim.integrator.s_points
            }
Rr = {"type" : "Relay",
         "value" : sim.relay.resp(),
         "t_points" : sim.integrator.t_points,
         "spatial_vec" : sim.integrator.s_points
        }

Ri = {"type" : "Interneuron",
         "value" : sim.interneuron.resp(),
         "t_points" : sim.integrator.t_points,
         "spatial_vec" : sim.integrator.s_points
}
Rc = {"type" : "Cortical",
            "value" : sim.cortical.resp(),
             "t_points" : sim.integrator.t_points,
             "spatial_vec" : sim.integrator.s_points
            }

data = [S, Rg, Rr, Ri,  Rc]
animate_imshow_plots(data, sim.integrator.dt,
colorbar = True, remove_axes = False,
save_animation = False, animation_name = "newTest")
