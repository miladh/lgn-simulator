#!/usr/bin/python
import os, sys
from sys import argv
from argparse import ArgumentParser
current_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(current_path, "..","..","tools"))
sys.path.append(lib_path)

import edog_runner

parser = ArgumentParser()
parser.add_argument("config_file", default=None)
args = parser.parse_args()

app_name = "nonlaggedXCells"

# Run edog:
config_file = args.config_file
run_id = edog_runner.run_edog(app_name, config_file)
