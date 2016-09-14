#!/usr/bin/python
import yaml
import os, os.path

def get_output_dir(record_label):
    from sumatra.projects import load_project
    data_path = os.path.abspath(load_project().data_store.root)
    output_dir = os.path.join(data_path, record_label)

    print "Results saved to this directory:\n", output_dir + "/*"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    return output_dir
