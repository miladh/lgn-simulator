from subprocess import call
from shutil import copyfile
import os, os.path
import yaml
import numpy as np

current_path = os.path.dirname(os.path.realpath(__file__))
copyfile(os.path.abspath(os.path.join(current_path,"phaseReversedPushPull.yaml")),
         os.path.abspath(os.path.join(current_path,"tmp.yaml")))
config_file = os.path.abspath(os.path.join(current_path,"tmp.yaml"))

with open(config_file, 'r') as stream:
    config_data = yaml.load(stream)

def modify_diameter(d):
    config_data["stimulus"]["maskSize"] = float(d)

def modify_kpg(kpg):
    config_data["stimulus"]["kId"] = int(kpg)

if __name__ == "__main__":
    spot_diameters = [0.07874015748031496, 0.31496062992125984 ]
    kpg = [3, 8]

    reason = "test if there is a temporal delay for low kpg"

    for k in kpg:
        modify_kpg(k)
        for d in spot_diameters:
            modify_diameter(d)

            with open(config_file, 'w') as stream:
                yaml.dump(config_data, stream)

            tag = "temporal_delay"

            call(["smt", "run", os.path.basename(config_file), "-i"+config_file, "-r "+ reason, "-t" +tag])

    os.remove(config_file)
