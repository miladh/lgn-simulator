import subprocess
import numpy as np

spot_diameters = np.linspace(0, 0.9, 1)

for d in spot_diameters:
    process = subprocess.Popen(["smt", "run", "spatialSummation.yaml", "maskSize="+str(d)],
                stdout=subprocess.PIPE)
    print process.communicate()[0]
