import subprocess
import numpy as np

spot_diameters = np.linspace(0, 0.9, 1)

for d in spot_diameters:
    process = subprocess.call(["smt", "run", "spatialSummation.yaml", "maskSize="+str(d)])
    print process.communicate()[0]
