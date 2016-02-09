from subprocess import call
import os
import os.path
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("project_dir", help='project directory')
parser.add_argument("root_data_path", help='root data path')


args = parser.parse_args()
root_data_path = args.root_data_path
project_dir = args.project_dir

if not os.path.exists(project_dir):
    print "project_dir not found: ", project_dir
    exit()

if not os.path.exists(root_data_path):
    print "root_data_path not found: ", root_data_path
    exit()


project_name = os.path.basename(project_dir)
data_path = os.path.join(root_data_path, project_name)
if not os.path.exists(data_path):
    os.makedirs(data_path)

current_path = os.path.dirname(os.path.realpath(__file__))
main_script = os.path.join(current_path, "run_app.py")


call("cd "+project_dir, shell=True)

if os.path.exists(".smt"):
    print "Sumatra project already initialized!"
else:
    print "initializing Sumatra project: ", project_name
    print "project directory: ", os.path.realpath(project_dir)
    print "data path: ", data_path

    call(["smt", "init", project_name, "-e", "python", "-m", main_script,
    "-l", "parameters", "-d", data_path])
