#!/usr/bin/python
import yaml
import subprocess
import os, os.path
from sumatra.projects import load_project

def run_simulator(config_file, record_label, run_id):

    current_path = os.path.dirname(os.path.realpath(__file__))
    app_name = load_project().name
    config_file = os.path.abspath(os.path.join(current_path, "../../apps", app_name, config_file))
    print "app_name: ", app_name
    print "config file: ", config_file

    with open(config_file, 'r') as stream:
        config_data = yaml.load(stream)
        # config_data["sumatra_label"] = str(config_data["sumatra_label"])
        # run_id = config_data["sumatra_label"]

    output_dir = os.path.join(os.path.abspath(load_project().data_store.root), record_label+ "_"+run_id)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        # os.makedirs(os.path.join(output_dir, record_label+ "_"+run_id))

    output_file = os.path.join(output_dir, record_label+ "_"+run_id + ".h5")
    print "output_file: ", output_file

    config_data["OutputManager"]["outputFilename"] = unicode(output_file)
    with open(config_file, 'w') as stream:
        yaml.safe_dump(config_data, stream, default_flow_style=False)

    build_path = os.path.abspath(os.path.join(current_path, "../../../build"))
    project_path = os.path.abspath(os.path.join(current_path, "../.."))

    print "Building in:\n", build_path

    if not os.path.exists(build_path):
        os.makedirs(build_path)
    subprocess.call(["qmake", project_path], cwd=build_path)
    subprocess.call(["make", "-j", "8"], cwd=build_path)

    app_path = os.path.join(build_path, "apps" , app_name)
    lib_path = os.path.join(build_path, "lib")

    env = dict(os.environ)
    env['LD_LIBRARY_PATH'] = lib_path

    run_argument = ["./lgnSimulator_"+ app_name, config_file, output_dir]
    print " ".join(run_argument)
    proc = subprocess.call(run_argument, cwd=app_path, env=env)

    print "Results saved to this directory:\n", output_dir + "/*"

    return run_id


# if __name__ == '__main__':
#     from argparse import ArgumentParser
#     parser = ArgumentParser()
#     parser.add_argument("config_file", help="app config file")
#     parser.add_argument("label", help="record label")
#     args = parser.parse_args()
#     config_file = args.config_file
#     label = args.label
#
#     print "Running lgn simulator..."
#     run_id = run_simulator(config_file, label)
