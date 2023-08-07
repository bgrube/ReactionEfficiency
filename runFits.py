#!/usr/bin/env python3


import os
import shutil
import subprocess


if __name__ == "__main__":

  # dataSample = "sig"
  # dataSample = "bggen"
  # dataSample = "data_030730"
  # dataSample = "data_041003"
  # dataSample = "data_042030"
  dataSample = "data_042550"
  studies = [
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_sigAllFixed_noBkg",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
      "pdfTypeBkg" : "\"\"",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_sigAllFudge_noBkg",
      "pdfTypeSig" : "Histogram",
      "pdfTypeBkg" : "\"\"",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_allFixed",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_sigSmear",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "shift scale",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_sigShift",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear scale",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_sigScale",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_sigFixSmear",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_sigFixShift",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "shift",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_sigFixScale",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "scale",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_sigAllFudge",
      "pdfTypeSig" : "Histogram",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift scale",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_bkgSmear",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "shift scale",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_bkgShift",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear scale",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_bkgScale",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear shift",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_bkgFixSmear",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "smear",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_bkgFixShift",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "shift",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_bkgFixScale",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
      "pdfTypeBkg" : "Histogram", "fixParsBkg" : "scale",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_bkgAllFudge",
      "pdfTypeSig" : "Histogram", "fixParsSig" : "smear shift scale",
      "pdfTypeBkg" : "Histogram",
    },
    {
      "fitDirectory" : f"BruFitOutput.{dataSample}_allFudge",
      "pdfTypeSig" : "Histogram",
      "pdfTypeBkg" : "Histogram",
    },
  ]


  for study in studies:
    fitDirectory = f"./fits/2018_01-ver02/noShowers/{study['fitDirectory']}"

    shutil.rmtree(fitDirectory, ignore_errors = True)
    os.makedirs(fitDirectory, exist_ok = True)
    print(f"Created directory '{fitDirectory}'")

    print(f"Starting fits ...")
    cmdLineOptions = [(f"--{option} {study[option]}" if option in study else "") for option in ("pdfTypeSig", "fixParsSig", "pdfTypeBkg", "fixParsBkg")]
    result = subprocess.run(f"./fitMissingMassSquared.py \"{fitDirectory}\" {' '.join(cmdLineOptions)} &> \"{fitDirectory}/fitMissingMassSquared.log\"", shell = True)
    if result.returncode != 0:
      raise RuntimeError(f"Fitting script failed with exit code '{result.returncode}'")

    subprocess.run(f"./cleanFitDir.sh \"{fitDirectory}\"", shell = True)
    subprocess.run(f"./cleanBruFitLogFile.py \"{fitDirectory}\"", shell = True)
    print("Plotting fit results...")
    subprocess.run(f"./plotFitResults.py \"{fitDirectory}\"  &> \"{fitDirectory}/plotFitResults.log\"", shell = True)
    print("Plotting efficiencies...")
    subprocess.run(f"./plotEfficiencies.py \"{fitDirectory}\"  &> \"{fitDirectory}/plotEfficiencies.log\"", shell = True)

    # pdfTypeBkg = "DoubleGaussian",
    # pdfTypeBkg = "DoubleGaussian_SameMean",
    # pdfTypeBkg = "SkewedGaussian_SkewNormal",
    # pdfTypeBkg = "SkewedGaussian_ExpMod",
    # pdfTypeBkg = "SkewedGaussian_Log",
