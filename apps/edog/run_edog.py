from subprocess import call
from shutil import copyfile
import os, os.path
import yaml
import numpy as np
import sys

current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../tools")))
import sumatra_tracking.run_simulator as st


def modify_num_spatial_points(ns):
    config_data["grid"]["ns"] = int(ns)

def modify_spatial_resolution(ds):
    config_data["grid"]["ds"] = float(ds)


def modify_phase(phase):
    config_data["stimulus"]["phase"] = float(phase)

def modify_spatial_freq(d):
    config_data["stimulus"]["kId"] = int(d)

def modify_diameter(d):
    config_data["stimulus"]["maskSize"] = float(d)



def modify_wrc(w):
    config_data["relay"]["Krc"]["w"] = float(w)

def modify_arc(a):
    config_data["relay"]["Krc"]["spatial"]["a"] = float(a)

def modify_brc(b):
    config_data["relay"]["Krc"]["spatial"]["b"] = float(b)

def modify_crc(c):
    config_data["relay"]["Krc"]["spatial"]["c"] = float(c)


def modify_arig(a):
    config_data["relay"]["Krig"]["spatial"]["a"] = float(a)

def modify_wrig(w):
    config_data["relay"]["Krig"]["w"] = float(w)



#read config file--------------------------------------------------------------
options = sys.argv[1:]
record_label = options[-1]
config_file = os.path.abspath(os.path.join(current_path, record_label+".yaml"))
copyfile(os.path.abspath(os.path.join(current_path,"edog.yaml")), config_file)

with open(config_file, 'r') as stream:
    config_data = yaml.load(stream)



#parameters---------------------------------------------------------------------
# widths = np.linspace(0, 3, 30)
# weights = np.linspace(0, 0.6, 3)
# diameters = np.linspace(0., 15, 250)
# w_rc = np.linspace(0, 0.9, 2)
# w_rc_c = np.linspace(0, 3, 4)
# spatial_freqs = range(0, 90)

weights = [0.8, 0.8, 0.5]
# weights_c = [0.5, 1.5, 2.0]
widths = [0.2, 1.0, 2.5]
# widths_b = [0.9, 2.5, 0.5]

# ds_vec = [0.0125, 0.025, 0.05, 0.1, 0.2]
ns_vec = [7, 8, 9, 10, 11]

#run simulator--------------------------------------------------------------------
counter= 0

modify_wrc(0.0)
modify_crc(0.0)
modify_arc(0.1)
modify_brc(0.9)
modify_wrig(-0.5)
modify_arig(0.3)

modify_spatial_resolution(0.05)
# modify_diameter(3.5)
# modify_spatial_freq(0)


for ns in ns_vec:
    modify_num_spatial_points(ns)
    for w, a in zip(weights, widths):
        modify_wrig(-w)
        modify_arig(a)
##########################################################
        with open(config_file, 'w') as stream:
            yaml.dump(config_data, stream)

        run_id = '{0:04}'.format(counter)
        st.run_simulator(config_file, record_label, run_id)
        counter+=1
os.remove(config_file)
##########################################################
