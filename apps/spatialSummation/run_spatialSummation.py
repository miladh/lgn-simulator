from subprocess import call
from shutil import copyfile
import os, os.path
import yaml
import numpy as np
import sys

current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../tools")))
import sumatra_tracking.run_simulator as st


def modify_spatial_freq(d):
    config_data["stimulus"]["kId"] = int(d)

def modify_diameter(d):
    config_data["stimulus"]["maskSize"] = float(d)

def modify_Kic(w):
    config_data["interneuron"]["Kic"]["w"] = float(w)

def modify_Krc(w):
    config_data["relay"]["Krc"]["w"]= float(w)

def modify_Kri(w):
    config_data["relay"]["Kri"]["w"]= float(-w)

#read config file-----------------------------------------------------------------------------
options = sys.argv[1:]
record_label = options[-1]
config_file = os.path.abspath(os.path.join(current_path, record_label+".yaml"))
copyfile(os.path.abspath(os.path.join(current_path,"spatialSummation.yaml")), config_file)

with open(config_file, 'r') as stream:
    config_data = yaml.load(stream)

#parameters-------------------------------------------------------------------------------------
spot_diameters = np.linspace(0., 6., 30)
# spatial_freqs = np.linspace(0.0, 120, 20)
weights = range(6)

#run simulator----------------------------------------------------------------------------------
counter= 0
for Kc in weights:
    modify_Krc(Kc)
    modify_Kic(Kc)
    for d in spot_diameters:
        modify_diameter(d)
        with open(config_file, 'w') as stream:
            yaml.dump(config_data, stream)

        run_id = '{0:04}'.format(counter)
        st.run_simulator(config_file, record_label, run_id)
        counter+=1

os.remove(config_file)
