#!/usr/bin/env python3


import os
import shutil
import subprocess


if __name__ == "__main__":

  bggenFileName = f"./pippippimpimpmiss_flatTree.MCbggen_2018_01-ver02.root.brufit"
  dataFileNames = [f"./pippippimpimpmiss_flatTree.RD_2018_01-ver02_041003.root.brufit"]
  dataSamples   = [{"dataFileName" : fileName, "dataLabel" : f"data_{fileName.split('.')[2].split('_')[-1]}"} for fileName in dataFileNames]
  print(f"!!! {dataSamples}")
  fits = [[
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_sigAllFixed_noBkg",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "\"\"",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_sigAllFudge_noBkg",
    #   "pdfTypeSig" : "Histogram",
    #   "pdfTypeBkg" : "\"\"",
    # },
    {
      "dataFileName" : dataSample['dataFileName'],
      "fitDirectory" : f"BruFitOutput.{dataSample['dataLabel']}_allFixed.new",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_sigSmear",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_sigShift",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_sigScale",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_sigFixSmear",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_sigFixShift",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "shift",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_sigFixScale",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_sigAllFudge",
    #   "pdfTypeSig" : "Histogram",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_bkgSmear",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "shift scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_bkgShift",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_bkgScale",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_bkgFixSmear",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_bkgFixShift",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "shift",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_bkgFixScale",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram", "fixParsBkg" : "scale",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_bkgAllFudge",
    #   "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
    #   "pdfTypeBkg" : "Histogram",
    # },
    # {
    #   "fitDirectory" : f"BruFitOutput.{dataSample}_allFudge",
    #   "pdfTypeSig" : "Histogram",
    #   "pdfTypeBkg" : "Histogram",
    # },
  ] for dataSample in dataSamples]
  print(f"!!! {fits}")
  # raise ValueError

  # pdfTypeBkg = "DoubleGaussian",
  # pdfTypeBkg = "DoubleGaussian_SameMean",
  # pdfTypeBkg = "SkewedGaussian_SkewNormal",
  # pdfTypeBkg = "SkewedGaussian_ExpMod",
  # pdfTypeBkg = "SkewedGaussian_Log",


  for dataSamples in fits:
    for study in dataSamples:
      fitDirectory = f"./fits/2018_01-ver02/noShowers/{study['fitDirectory']}"
      # prepare directories
      shutil.rmtree(fitDirectory, ignore_errors = True)
      os.makedirs(fitDirectory, exist_ok = True)
      print(f"Created directory '{fitDirectory}'")
      # run fits
      print(f"Starting fits ...")
      cmdLineOptions = [(f"--{option} {study[option]}" if option in study else "") for option in ("pdfTypeSig", "fixParsSig", "pdfTypeBkg", "fixParsBkg")]
      result = subprocess.run(
        f"./fitMissingMassSquared.py \"{fitDirectory}\" {study['dataFileName']} {bggenFileName} {' '.join(cmdLineOptions)} &> \"{fitDirectory}/fitMissingMassSquared.log\"", shell = True)
      if result.returncode != 0:
        raise RuntimeError(f"Fitting script failed with exit code '{result.returncode}'")
      # postprocess fit results
      subprocess.run(f"./cleanFitDir.sh \"{fitDirectory}\"", shell = True)
      subprocess.run(f"./cleanBruFitLogFile.py \"{fitDirectory}\"", shell = True)
      print("Plotting fit results...")
      subprocess.run(f"./plotFitResults.py \"{fitDirectory}\"  &> \"{fitDirectory}/plotFitResults.log\"", shell = True)
      print("Plotting efficiencies...")
      subprocess.run(f"./plotEfficiencies.py \"{fitDirectory}\"  &> \"{fitDirectory}/plotEfficiencies.log\"", shell = True)
