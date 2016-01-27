from subprocess import call
import numpy as np

spot_diameters = np.linspace(0, 0.8, 20)

for d in spot_diameters:
    call(["smt", "run", "nonlaggedXCellsConfig.yaml", "maskSize="+str(d)])
