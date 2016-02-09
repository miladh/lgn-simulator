#!/usr/bin/python
import os, sys
from sys import argv
from argparse import ArgumentParser
current_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(current_path, "..","..","tools"))
sys.path.append(lib_path)

import lgnSimulator_runner

parser = ArgumentParser()
parser.add_argument("config_file", default=None)
args = parser.parse_args()

app_name = os.path.basename(current_path)

# Run edog:
config_file = args.config_file
run_id = lgnSimulator_runner.run_lgnSimulator(app_name, config_file)
