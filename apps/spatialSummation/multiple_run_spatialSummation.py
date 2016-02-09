from subprocess import call
import numpy as np

call(["smt", "run", "spatialSummation.yaml", "spatialSummation",
"maskSize="+str(0.3)])
