from subprocess import call
import numpy as np

call(["smt", "run", "spatialSummation.yaml",  "maskSize="+str(0.1)])
