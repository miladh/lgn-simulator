from subprocess import call
from shutil import copyfile
import os, os.path
import yaml
import numpy as np
import sys

current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../tools")))
import sumatra_tracking.run_simulator as st


def run_simulation(attrs):
    counter= 0
    print attrs
    for key, attr in attrs.iteritems():
        print attr, key
        for value in attr["values"]:
            attr["func"](value)
            print value
            with open(config_file, 'w') as stream:
                yaml.dump(config_data, stream)

            run_id = '{0:04}'.format(counter)
            st.run_simulator(config_file, record_label, run_id)
            counter+=1

    os.remove(config_file)


def modify_phase(phase):
    config_data["stimulus"]["phase"] = float(phase)

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

def modify_arc(a):
    config_data["relay"]["Krc"]["spatial"]["a"] = float(a)



#read config file-----------------------------------------------------------------------------
options = sys.argv[1:]
record_label = options[-1]
config_file = os.path.abspath(os.path.join(current_path, record_label+".yaml"))
copyfile(os.path.abspath(os.path.join(current_path,"spatialSummation.yaml")), config_file)

with open(config_file, 'r') as stream:
    config_data = yaml.load(stream)


# #parameters-------------------------------------------------------------------------------------
spot_diameters = np.linspace(0., 15, 1)
spatial_freqs = 4
weights = np.linspace(0, 1.0, 6)
phase = np.linspace(-180, 180, 10)
widths = np.array([0.25, 0.5, 1.0, 2.0])*0.1

attrs = {
    "a_rc" : {"func": modify_arc, "values": widths}
}

#run simulator----------------------------------------------------------------------------------
run_simulation(attrs)
