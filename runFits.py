#!/usr/bin/env python3


import functools
import glob
import os
import shutil
import subprocess
from typing import (
  Dict,
  List,
)


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


if __name__ == "__main__":

  dataPeriods = (
    "2017_01-ver03",
    "2018_01-ver02",
    "2018_08-ver02",
    "2019_11-ver01",
  )

  dataSamples: List[Dict[str, str]] = []
  for dataPeriod in dataPeriods:
    dataSamples += [
      { # bggen MC
        **dict.fromkeys(["dataFileName", "bggenFileName"],
                        f"./data/MCbggen/{dataPeriod}/pippippimpimpmiss_flatTree.MCbggen_{dataPeriod}.root.brufit"),  # dataFileName and bggenFileName have identical value
        "dataPeriod" : dataPeriod,
        "dataLabel"  : f"bggen_{dataPeriod}",
      },
      { # real data
        "dataFileName"  : f"./data/RD/{dataPeriod}/pippippimpimpmiss_flatTree.RD_{dataPeriod}.root.brufit",
        "bggenFileName" : f"./data/MCbggen/{dataPeriod}/pippippimpimpmiss_flatTree.MCbggen_{dataPeriod}.root.brufit",
        "dataPeriod"    : dataPeriod,
        "dataLabel"     : f"data_{dataPeriod}",
      },
    ]

  fits: List[List[Dict[str, str]]] = [[  # list of fits for each data sample
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigAllFixed_noBkg",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "\"\"",  #TODO is it correct to define a string '""'? why not just empty string?
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigAllFudge_noBkg",
    #   "pdfTypeSig" : "Histogram",
    #   "pdfTypeBkg" : "\"\"",
    # },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_allFixed",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigSmear",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigShift",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigScale",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigFixSmear",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigFixShift",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "shift",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigFixScale",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigAllFudge",
    #   "pdfTypeSig" : "Histogram",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgSmear",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgShift",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgScale",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgFixSmear",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgFixShift",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "shift",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgFixScale",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgAllFudge",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_allFudge",
    #   "pdfTypeSig" : "Histogram",
    #   "pdfTypeBkg" : "Histogram",
    # },
  ] for dataSample in dataSamples]
  # add data period and input files (same for all fits of a given data sample)
  for index, dataSample in enumerate(dataSamples):
    for fit in fits[index]:
      fit.update({"dataPeriod"    : dataSample["dataPeriod"]})
      fit.update({"dataFileName"  : dataSample["dataFileName"]})
      fit.update({"bggenFileName" : dataSample["bggenFileName"]})

  # pdfTypeBkg = "DoubleGaussian",
  # pdfTypeBkg = "DoubleGaussian_SameMean",
  # pdfTypeBkg = "SkewedGaussian_SkewNormal",
  # pdfTypeBkg = "SkewedGaussian_ExpMod",
  # pdfTypeBkg = "SkewedGaussian_Log",

  for fitsForDataSample in fits:
    for fit in fitsForDataSample:
      fitDirectory = f"./fits/{fit['dataPeriod']}/noShowers/{fit['fitDirectory']}"
      # prepare directories
      shutil.rmtree(fitDirectory, ignore_errors = True)
      os.makedirs(fitDirectory, exist_ok = True)
      print(f"Created directory '{fitDirectory}'")
      # run fits
      #TODO call python functions directly instead of making the detour via the command-line interface
      print(f"Starting fits ...")
      cmdLineOptions = [(f"--{option} {fit[option]}" if option in fit else "") for option in ("pdfTypeSig", "fixParsSig", "pdfTypeBkg", "fixParsBkg")]
      result = subprocess.run(
        f"./fitMissingMassSquared.py \"{fitDirectory}\" {fit['dataFileName']} {fit['bggenFileName']} {' '.join(cmdLineOptions)} &> \"{fitDirectory}/fitMissingMassSquared.log\"", shell = True)
      if result.returncode != 0:
        raise RuntimeError(f"Fitting script failed with exit code '{result.returncode}'")
      # postprocess fit results
      subprocess.run(f"./cleanFitDir.sh \"{fitDirectory}\"", shell = True)
      subprocess.run(f"./cleanBruFitLogFile.py \"{fitDirectory}\"", shell = True)
      print("Plotting fit results...")
      subprocess.run(f"./plotFitResults.py \"{fitDirectory}\"  &> \"{fitDirectory}/plotFitResults.log\"", shell = True)
      print("Plotting efficiencies...")
      subprocess.run(f"./plotEfficiencies.py \"{fitDirectory}\"  &> \"{fitDirectory}/plotEfficiencies.log\"", shell = True)
