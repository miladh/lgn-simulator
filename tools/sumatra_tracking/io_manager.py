#!/usr/bin/python
import yaml
import os, os.path

def get_output_dir(config_file):
    from sumatra.projects import load_project
    with open(config_file, 'r') as stream:
        config_data = yaml.load(stream)


    run_id = config_data["sumatra_label"]
    data_path = os.path.abspath(load_project().data_store.root)
    output_dir = os.path.join(data_path, run_id)

    print "Results saved to this directory:\n", output_dir + "/*"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    return output_dir


def get_simulations(config_file):
    import h5py
    from analysis.Simulation import Simulation

    with open(config_file, 'r') as stream:
        config_data = yaml.load(stream)
        data_path = config_data["data_path"]
        data_ids  = config_data["simulation_ids"]

    print "Data path: ", data_path
    print "Reading simulation_ids:\n", data_ids

    sims=[]
    for data_id in data_ids:
        data_file = os.path.join(data_path, data_id, data_id +".h5")
        f = h5py.File(data_file, "r")
        sims.append(Simulation(f))

    return sims
