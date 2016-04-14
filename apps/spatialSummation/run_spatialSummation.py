from subprocess import call
from shutil import copyfile
import os, os.path
import yaml
import numpy as np

current_path = os.path.dirname(os.path.realpath(__file__))
copyfile(os.path.abspath(os.path.join(current_path,"spatialSummation.yaml")),
         os.path.abspath(os.path.join(current_path,"tmp.yaml")))
config_file = os.path.abspath(os.path.join(current_path,"tmp.yaml"))

with open(config_file, 'r') as stream:
    config_data = yaml.load(stream)

def modify_diameter(d):
    config_data["stimulus"]["maskSize"] = str(d)

def modify_inhibition_weight(w):
    config_data["interneuron"]["Kic"]["spatial"]["A"] = str(w)

if __name__ == "__main__":
    spot_diameters = np.linspace(0, 0.9, 1)
    weights = np.linspace(0.1, 2, 1)

    reason = "Test the effect inhibition weight on area summation curves"

    for w in weights:
        modify_inhibition_weight(w)
        for d in spot_diameters:
            modify_diameter(d)

            with open(config_file, 'w') as stream:
                yaml.dump(config_data, stream)

            tag = "d="+str(d)+ ", w_inhib="+str(2)

            call(["smt", "run", os.path.basename(config_file), "-r "+ reason, "-t" +tag])

    # os.remove(config_file)
