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
    config_data["stimulus"]["maskSize"] = float(d)

def modify_Kic(w):
    config_data["interneuron"]["Kic"]["spatial"]["A"] = float(w)

def modify_Krc(w):
    config_data["relay"]["Krc"]["spatial"]["A"] = float(w)

def modify_Kri(w):
    config_data["relay"]["Kri"]["spatial"]["weight"]= float(-w)

if __name__ == "__main__":
    spot_diameters = np.linspace(0, 0.9, 4)
    weights = np.linspace(0.1, 2.0, 4)

    reason = "Test run with new setup"

    for w in weights:
        modify_Kic(w)
        modify_Krc(w)

        # modify_Kri(w)
        for d in spot_diameters:
            modify_diameter(d)

            with open(config_file, 'w') as stream:
                yaml.dump(config_data, stream)

            tag = "test_run_3"

            call(["smt", "run", os.path.basename(config_file), "-i"+config_file, "-r "+ reason, "-t" +tag])

    os.remove(config_file)
