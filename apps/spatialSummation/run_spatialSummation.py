from subprocess import call
from shutil import copyfile
import os, os.path
import yaml
import numpy as np
import sys

current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../tools")))
import sumatra_tracking.run_simulator as st


def modify_diameter(d):
    config_data["stimulus"]["maskSize"] = float(d)

def modify_Kic(w):
    config_data["interneuron"]["Kic"]["spatial"]["A"] = float(w)

def modify_Krc(w):
    config_data["relay"]["Krc"]["spatial"]["A"] = float(w)

def modify_Kri(w):
    config_data["relay"]["Kri"]["spatial"]["weight"]= float(-w)

#read config file
options = sys.argv[1:]
record_label = options[-1]
copyfile(os.path.abspath(os.path.join(current_path,"spatialSummation.yaml")),
         os.path.abspath(os.path.join(current_path, record_label+".yaml")))
config_file = os.path.abspath(os.path.join(current_path, record_label+".yaml"))

with open(config_file, 'r') as stream:
    config_data = yaml.load(stream)

#parameters
spot_diameters = np.linspace(0, 0.9, 2)
weights = np.linspace(0.0, 2.0, 10)



counter= 0
# for w in weights:
#     modify_Kic(w)
#     modify_Krc(w)

    # modify_Kri(w)

for d in spot_diameters:
    modify_diameter(d)

    with open(config_file, 'w') as stream:
        yaml.dump(config_data, stream)

    run_id = '{0:04}'.format(counter)
    tmp_config_filename = record_label+"_"+run_id +".yaml"

    copyfile(config_file, os.path.abspath(os.path.join(current_path, tmp_config_filename)))
    tmp_config_file = os.path.abspath(os.path.join(current_path, tmp_config_filename))

    st.run_simulator(tmp_config_file, record_label, run_id)
    os.remove(tmp_config_file)
    counter+=1


os.remove(config_file)
