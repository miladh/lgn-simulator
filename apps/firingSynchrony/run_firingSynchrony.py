from subprocess import call
from shutil import copyfile
import os, os.path
import yaml
import numpy as np
import sys

current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../tools")))
import sumatra_tracking.run_simulator as st

def modify_phase(phase):
    config_data["stimulus"]["phase"] = float(phase)

def modify_spatial_freq(d):
    config_data["stimulus"]["kId"] = int(d)

def modify_temp_freq(d):
    config_data["stimulus"]["wId"] = int(d)

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


def modify_tau_rc(t):
    config_data["relay"]["Krc"]["temporal"]["tau"] = float(t)

def modify_delay_rc(t):
    config_data["relay"]["Krc"]["temporal"]["delay"] = float(t)



def modify_arig(a):
    config_data["relay"]["Krig"]["spatial"]["a"] = float(a)

def modify_wrig(w):
    config_data["relay"]["Krig"]["w"] = float(w)

def modify_tau_rig(t):
    config_data["relay"]["Krig"]["temporal"]["tau"] = float(t)

def modify_delay_rig(t):
    config_data["relay"]["Krig"]["temporal"]["delay"] = float(t)

#read config file--------------------------------------------------------------
options = sys.argv[1:]
record_label = options[-1]
config_file = os.path.abspath(os.path.join(current_path, record_label+".yaml"))
copyfile(os.path.abspath(os.path.join(current_path,"firingSynchrony.yaml")), config_file)

with open(config_file, 'r') as stream:
    config_data = yaml.load(stream)


#parameters---------------------------------------------------------------------
weights = np.linspace(0, 0.9, 4)

tau_vec = np.linspace(1,50,15)
delay_vec = np.linspace(4,32,16)

modify_arc(0.1)
modify_brc(0.9)
modify_crc(2.0)
modify_wrc(0.5)

modify_wrig(-0.0)
modify_arig(0.3)

# modify_tau_rc(1)
# modify_delay_rc(0)

modify_tau_rig(1)
modify_delay_rig(0)


#run simulator--------------------------------------------------------------------
counter= 0
for tau in tau_vec:
    modify_tau_rc(tau)
    for delay in delay_vec:
        modify_delay_rc(delay)
##########################################################
        with open(config_file, 'w') as stream:
            yaml.dump(config_data, stream)

        run_id = '{0:04}'.format(counter)
        st.run_simulator(config_file, record_label, run_id)
        counter+=1
os.remove(config_file)
##########################################################
