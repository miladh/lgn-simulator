#!/usr/bin/python
import yaml
import os, os.path

def get_output_dir(config_file, run_id):
    from sumatra.projects import load_project
    with open(config_file, 'r') as stream:
        config_data = yaml.load(stream)

    data_path = os.path.abspath(load_project().data_store.root)
    output_dir = os.path.join(data_path, run_id)

    print "Results saved to this directory:\n", output_dir + "/*"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    return output_dir
