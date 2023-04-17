#!/usr/bin/env python3


import os
import shutil
import subprocess


if __name__ == "__main__":

  fitDirectory = "./fits/noShowers/BruFitOutput.data_sigAllFudge_bkgSkewedGaussian_ExpMod"

  shutil.rmtree(fitDirectory, ignore_errors = True)
  os.makedirs(fitDirectory, exist_ok = True)
  print(f"Created directory '{fitDirectory}'")

  print(f"Starting fits ...")
  result = subprocess.run(f"./fitMissingMassSquared.py \"{fitDirectory}\" &> \"{fitDirectory}/fitMissingMassSquared.log\"", shell = True)
  if result.returncode != 0:
    raise RuntimeError(f"Fitting script failed with exit code '{result.returncode}'")

  subprocess.run(f"./cleanFitDir.sh \"{fitDirectory}\"", shell = True)
  subprocess.run(f"./cleanBruFitLogFile.py \"{fitDirectory}\"", shell = True)
  print("Plotting fit results...")
  subprocess.run(f"./plotFitResults.py \"{fitDirectory}\"  &> \"{fitDirectory}/plotFitResults.log\"", shell = True)
  print("Plotting efficiencies...")
  subprocess.run(f"./plotEfficiencies.py \"{fitDirectory}\"  &> \"{fitDirectory}/plotEfficiencies.log\"", shell = True)
