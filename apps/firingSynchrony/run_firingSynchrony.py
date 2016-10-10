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


def modify_arc(a):
    config_data["relay"]["Krc"]["spatial"]["a"] = float(a)

def modify_wrc(w):
    config_data["relay"]["Krc"]["w"] = float(w)

def modify_tau_rc(t):
    config_data["relay"]["Krc"]["temporal"]["tau"] = float(t)

def modify_delay_rc(t):
    config_data["relay"]["Krc"]["temporal"]["delay"] = float(t)




def modify_wcr(w):
    config_data["cortical"]["Kcr"]["w"] = float(w)



def modify_aic(a):
    config_data["interneuron"]["Kic"]["spatial"]["a"] = float(a)

def modify_wic(w):
    config_data["interneuron"]["Kic"]["w"] = float(w)


def modify_aig(a):
    config_data["interneuron"]["Kig"]["spatial"]["a"] = float(a)

def modify_wig(w):
    config_data["interneuron"]["Kig"]["w"] = float(w)


def modify_ari(a):
    config_data["relay"]["Kri"]["spatial"]["a"] = float(a)

def modify_wri(w):
    config_data["relay"]["Kri"]["w"] = float(w)

def modify_tau_ri(t):
    config_data["relay"]["Kri"]["temporal"]["tau"] = float(t)

def modify_delay_ri(t):
    config_data["relay"]["Kri"]["temporal"]["delay"] = float(t)

#read config file--------------------------------------------------------------
options = sys.argv[1:]
record_label = options[-1]
config_file = os.path.abspath(os.path.join(current_path, record_label+".yaml"))
copyfile(os.path.abspath(os.path.join(current_path,"firingSynchrony.yaml")), config_file)

with open(config_file, 'r') as stream:
    config_data = yaml.load(stream)


#parameters---------------------------------------------------------------------
w_ic = [0, 0.8, 2.4, 4.0]
w_rc = [0, 0.18, 0.54, 0.9]
spatial_freqs = range(0, 60)

tau_ri = np.linspace(1,50,15)
delay_ri = np.linspace(2,32,16)

modify_wic(0)
modify_wrc(0.5)
modify_wcr(1)
modify_wri(0)
modify_wig(0)

#run simulator--------------------------------------------------------------------
counter= 0
for tau in tau_ri:
    modify_tau_rc(tau)
    for delay in delay_ri:
        modify_delay_rc(delay)
##########################################################
        with open(config_file, 'w') as stream:
            yaml.dump(config_data, stream)

        run_id = '{0:04}'.format(counter)
        st.run_simulator(config_file, record_label, run_id)
        counter+=1
os.remove(config_file)
##########################################################
