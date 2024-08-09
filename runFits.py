#!/usr/bin/env python3


from __future__ import annotations

import functools
import os
import shutil
import subprocess
from typing import Any
from wurlitzer import pipes, STDOUT

import ROOT

from fitMissingMassSquared import fitMissingMassSquared


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


if __name__ == "__main__":

  fitRootDir = "./fits"
  # fitRootDir = "./fits.pionComparison"
  # fitRootDir = "./fits.pionComparison.R6.28"
  dataPeriods = (
    # "2017_01-ver03",
    # "2018_01-ver02",
    # "2018_08-ver02",
    # "2019_11-ver01",
    "2017_01-ver03_goodToF",
    # "2018_01-ver02_goodToF",
    # "2018_08-ver02_goodToF",
  )
  useTotal = False

  dataSamples: list[dict[str, str]] = []
  for dataPeriod in dataPeriods:
    dataSamples += [
      # { # bggen MC
      #   **dict.fromkeys(["dataFileName", "bggenFileName"],
      #                   f"./data/MCbggen/{dataPeriod}/pippippimpimpmiss_flatTree.MCbggen_{dataPeriod}.root.brufit"),  # "dataFileName" and "bggenFileName" are set to identical values
      #   "dataPeriod" : dataPeriod,
      #   "dataLabel"  : f"bggen_{dataPeriod}",
      # },
      { # real data
        "dataFileName"  : f"./data/RD/{dataPeriod}/pippippimpimpmiss_flatTree.RD_{dataPeriod}.root.brufit",
        "bggenFileName" : f"./data/MCbggen/{dataPeriod}/pippippimpimpmiss_flatTree.MCbggen_{dataPeriod}.root.brufit",
        "dataPeriod"    : dataPeriod,
        "dataLabel"     : f"data_{dataPeriod}",
      },
    ]

  fits: list[list[dict[str, Any]]] = [[  # list of fits for each data sample
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigAllFixed_noBkg",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #     "pdfTypeBkg" : "\"\"",  #TODO is it correct to define a string '""'? why not just empty string?
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigAllFudge_noBkg",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram",
    #     "pdfTypeBkg" : "\"\"",
    #   },
    # },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_allFixed",
      "kwargs" : {
        "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
        "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
      },
    },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigSmear",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "shift scale",
    #     "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigShift",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "smear scale",
    #     "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigScale",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift",
    #     "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigFixSmear",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "smear",
    #     "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigFixShift",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "shift",
    #     "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigFixScale",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "scale",
    #     "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_sigAllFudge",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram",
    #     "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgSmear",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #     "pdfTypeBkg" : "Histogram", "fixParsBkg" : "shift scale",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgShift",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #     "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear scale",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgScale",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #     "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgFixSmear",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #     "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgFixShift",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #     "pdfTypeBkg" : "Histogram", "fixParsBkg" : "shift",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgFixScale",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #     "pdfTypeBkg" : "Histogram", "fixParsBkg" : "scale",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_bkgAllFudge",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #     "pdfTypeBkg" : "Histogram",
    #   },
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_allFudge",
    #   "kwargs" : {
    #     "pdfTypeSig" : "Histogram",
    #     "pdfTypeBkg" : "Histogram",
    #   },
    # },
  ] for dataSample in dataSamples]
  # add data period and input files (same for all fits of a given data sample)
  for index, dataSample in enumerate(dataSamples):
    for fit in fits[index]:
      fit.update({"dataPeriod" : dataSample["dataPeriod"]})
      fit["kwargs"]["dataFileName"]  = dataSample["dataFileName"]
      fit["kwargs"]["bggenFileName"] = dataSample["bggenFileName"]

  ROOT.gBenchmark.Start("Total processing time")
  for fitsForDataSample in fits:
    for fit in fitsForDataSample:
      fitDirectory = f"{fitRootDir}/{fit['dataPeriod']}/noShowers/{fit['fitDirectory']}"
      # recreate directories if already existing
      shutil.rmtree(fitDirectory, ignore_errors = True)
      os.makedirs(fitDirectory, exist_ok = True)
      print(f"Created directory '{fitDirectory}'")
      # run fits
      print(f"Starting fits ...")
      with open(f"{fitDirectory}/fitMissingMassSquared.log", "w") as logFile, pipes(logFile, stderr = STDOUT):  # write separate log file for each fit
        fitMissingMassSquared(
          outputDirName = fitDirectory,
          **fit["kwargs"],
        )
      # postprocess fit results
      print("Plotting fit results...")
      #TODO call python functions directly instead of making the detour via the command-line interface
      subprocess.run(f"./plotFitResults.py \"{fitDirectory}\"  &> \"{fitDirectory}/plotFitResults.log\"", shell = True)
      print("Plotting efficiencies...")
      subprocess.run(f"./plotEfficiencies.py {'--useTotal' if useTotal else ''} \"{fitDirectory}\"  &> \"{fitDirectory}/plotEfficiencies.log\"", shell = True)
  ROOT.gBenchmark.Show("Total processing time")
