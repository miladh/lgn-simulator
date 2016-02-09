from subprocess import call
import numpy as np

app_name = "spatialSummation"
spot_diameters = np.linspace(0, 0.9, 1)

for d in spot_diameters:
    call(["smt", "run", "spatialSummation.yaml", app_name,
    "maskSize="+str(d)])
