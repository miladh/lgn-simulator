#!/usr/bin/python
import yaml
import os, os.path
import h5py
from sumatra.projects import load_project
import Simulation

def get_simulation_environment(config_file, record):
    run_id = "tmp"
    with open(config_file, 'r') as stream:
        config_data = yaml.load(stream)
        data_ids = config_data["simulation_ids"]
        if record:
            run_id = config_data["sumatra_label"]

    data_path = os.path.abspath(load_project().data_store.root)
    output_dir = os.path.join(data_path, run_id)

    print "simulation_ids: ", data_ids
    print "Results saved to this directory:\n", output_dir + "/*"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read data:--------------------------------------------------------------------
    sims = []
    for data_id in data_ids:
        data_file = os.path.join(data_path, data_id, data_id +".h5")
        f = h5py.File(data_file, "r")
        sims.append(Simulation.Simulation(f))

    return sims, output_dir
