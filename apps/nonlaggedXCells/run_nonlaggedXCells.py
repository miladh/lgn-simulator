from subprocess import call
import numpy as np

spot_diameters = np.linspace(0, 0.9, 30)

for d in spot_diameters:
    call(["smt", "run", "nonlaggedXCells.yaml", "maskSize="+str(d)])
