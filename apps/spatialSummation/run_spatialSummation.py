from subprocess import call
import numpy as np

spot_diameters = np.linspace(0, 0.9, 30)
weights = np.linspace(0.1, 2, 5)

for w in weights:
    for d in spot_diameters:
        call(["smt", "run", "spatialSummation.yaml", "-t w= " + str(w),
        "inhibitionWeight="+str(w), "maskSize="+str(d)])
