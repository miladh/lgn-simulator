#!/usr/bin/python
import os, sys
from sys import argv
from argparse import ArgumentParser
current_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(current_path, "..","..","tools"))
sys.path.append(lib_path)

import Edog_runner as edog_runner

parser = ArgumentParser()
parser.add_argument("--id", nargs='?')
parser.add_argument("config_file")
args = parser.parse_args()

config_file = args.config_file
run_id = "spatialSummation"
app_name = "spatialSummation"

edog_runner.run_edog(args.id, app_name, config_file)

#Plotting
