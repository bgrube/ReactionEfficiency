#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import itertools
# import os

from uncertainties import ufloat
from uncertainties import umath

import ROOT

import makePlots  # defines helper functions to generate histograms from data trees

makePlots.setupPlotStyle()


YIELD_PAR_NAMES = {
  "Signal"     : "Yld_SigPdf",
  "Background" : "Yld_BgPdf"
}


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable

  # dataset           = "030730",
  dataset           = "bggen_2017_01-ver03"
  dataFileName      = f"../ReactionEfficiency/pippippimpimpmiss_flatTree.{dataset}.root.brufit"
  outputDirName     = "BruFitOutput"
  kinematicBinnings = [  # list of tuples [ (variable, nmb of bins, min value, max value) ]
    ("BeamEnergy", 9, 3.0, 12.0)
  ]

  dataSets = ["Total", "Found", "Missing"]
  dataSet = "Total"

  binningFileName = f"{outputDirName}/{dataSet}/DataBinsConfig.root"
  print(f"Loading binning from file '{binningFileName}'")
  bins = ROOT.Bins("HSBins", binningFileName)
  print("Found binning:")
  bins.PrintAxis()
  binFileNames = [str(fileName) for fileName in bins.GetFileNames()]
  # verify that all bin files exist
  for binFileName in binFileNames:
    print(binFileName)
  axes = bins.GetVarAxis()
  # for axis in axes:
  #   print(axis.GetName(), axis.GetNbins())
  #   for axisBinIndex in range(1, axis.GetNbins() + 1):
  #     axisBinCenter = axis.GetBinCenter(axisBinIndex)
  #     print(axis.GetName(), axisBinCenter)

  yields = []  # list of dictionaries with all yields [ { <bin var> : <bin center>, ..., "Signal" : <yield>, "Background" : <yield> } ]
  axisBinIndexRanges = tuple(range(1, axis.GetNbins() + 1) for axis in axes)
  for axisBinIndices in itertools.product(*axisBinIndexRanges):
    # for axisIndex, axisBinIndex in enumerate(axisBinIndices):
    #   print(axes[axisIndex].GetName(), axisBinIndex, axes[axisIndex].GetBinCenter(axisBinIndex), "   ", end = "")
    # print()
    # axisNames      = tuple(axes[axisIndex].GetName()                  for axisIndex, _            in enumerate(axisBinIndices))
    axisBinCenters    = tuple(axes[axisIndex].GetBinCenter(axisBinIndex) for axisIndex, axisBinIndex in enumerate(axisBinIndices))
    yieldsInBin       = {axes[axisIndex].GetName() : axisBinCenter for axisIndex, axisBinCenter in enumerate(axisBinCenters)}  # store bin info
    binIndex          = bins.FindBin(*axisBinCenters)  #!Note! the unpacking works only up to 6 binning dimensions
    # binName           = bins.GetBinName(binIndex)
    fitResultFileName = binFileNames[binIndex].replace("TreeData.root", "ResultsHSMinuit2.root")
    print(f"Reading fit result object 'MinuitResult' from file '{fitResultFileName}'")
    # print(axisNames, axisBinIndices, axisBinCenters, binIndex, binName, fitResultFileName)
    fitResultFile     = ROOT.TFile.Open(fitResultFileName, "READ")
    fitResult         = fitResultFile.Get("MinuitResult")
    fitPars           = fitResult.floatParsFinal()
    # for fitParIndex in range(fitPars.getSize()):
    #   fitPar = fitPars[fitParIndex]
    #   print(fitPar, "    ", end = "")
    # print()
    for yieldType, yieldParName in YIELD_PAR_NAMES.items():
      yieldParIndex = fitPars.index(yieldParName)
      if yieldParIndex < 0:
        raise IndexError(f"Cannot find parameter '{yieldParName}' in fit parameters {fitPars}")
      yieldPar = fitPars[yieldParIndex]
      yieldsInBin[yieldType] = ufloat(yieldPar.getVal(), yieldPar.getError())  # store yields
      # print(f"{yieldType} yield is: {yieldsInBin[yieldType]}")
    print(f"Read yields: {yieldsInBin}")
    yields.append(yieldsInBin)
  print(yields)
