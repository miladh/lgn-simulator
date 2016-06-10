from subprocess import call
from shutil import copyfile
import os, os.path
import yaml
import numpy as np
import sys

current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(current_path,"../../tools")))

import sumatra_tracking.run_app as st


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
    spot_diameters = np.linspace(0, 0.9, 2)
    weights = np.linspace(0.0, 2.0, 10)

    options = sys.argv[1:]
    record_label = options[-1]

    # for w in weights:
    #     modify_Kic(w)
    #     modify_Krc(w)

        # modify_Kri(w)

    counter= 0
    for d in spot_diameters:
        modify_diameter(d)

        with open(config_file, 'w') as stream:
            yaml.dump(config_data, stream)

        run_id = '{0:03}'.format(counter)
        st.run_simulator(os.path.basename(config_file), record_label, run_id)
        counter+=1
    os.remove(config_file)


# if __name__ == "__main__":
#     spot_diameters = np.linspace(0, 0.9, 10)
#     weights = np.linspace(0.0, 2.0, 10)
#
#     reason = "Kc effect, with different width for I and R"
#
#     for w in weights:
#         modify_Kic(w)
#         modify_Krc(w)
#
#         # modify_Kri(w)
#         for d in spot_diameters:
#             modify_diameter(d)
#
#             with open(config_file, 'w') as stream:
#                 yaml.dump(config_data, stream)
#
#             tag = "Kc_different_width_pg"
#
#             call(["smt", "run", os.path.basename(config_file), "-i"+config_file, "-r "+ reason, "-t" +tag])
#
#     os.remove(config_file)
