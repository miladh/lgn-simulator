from subprocess import call
from shutil import copyfile
import os, os.path
import yaml
import numpy as np

current_path = os.path.dirname(os.path.realpath(__file__))
copyfile(os.path.abspath(os.path.join(current_path,"stimuliAnalysis.yaml")),
         os.path.abspath(os.path.join(current_path,"tmp.yaml")))
config_file = os.path.abspath(os.path.join(current_path,"tmp.yaml"))

with open(config_file, 'r') as stream:
    config_data = yaml.load(stream)

def modify_diameter(d):
    config_data["stimulus"]["maskSize"] = float(d)


if __name__ == "__main__":
    spot_diameters = np.linspace(0, 0.9, 10)

    reason = "exploring FT of patch stim"

    for d in spot_diameters:
        modify_diameter(d)

        with open(config_file, 'w') as stream:
            yaml.dump(config_data, stream)

        tag = "stim_different_d"

        call(["smt", "run", os.path.basename(config_file), "-i"+config_file, "-r "+ reason, "-t" +tag])

os.remove(config_file)
