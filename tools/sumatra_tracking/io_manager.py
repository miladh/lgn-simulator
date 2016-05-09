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
