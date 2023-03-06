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


# returns
#     tuple with binning variables (<binning var>, ... )
#     list of dict with yields [ { <binning var> : <bin center>, ..., "Signal" : <yield>, "Background" : <yield> }, ... ]
def readYieldsFromFitResults(fitResultDirName):
  # get binning
  binningFileName = f"{fitResultDirName}/DataBinsConfig.root"
  print(f"Loading binning from file '{binningFileName}'")
  bins = ROOT.Bins("HSBins", binningFileName)
  print("Found binning:")
  bins.PrintAxis()
  axes         = bins.GetVarAxis()
  binVarNames  = tuple(axis.GetName() for axis in axes)
  binFileNames = [str(fileName) for fileName in bins.GetFileNames()]
  # read yields from files with fit results
  yields = []
  axisBinIndexRanges = tuple(range(1, axis.GetNbins() + 1) for axis in axes)
  for axisBinIndices in itertools.product(*axisBinIndexRanges):  # loop over all tuples of bin indices for the axes
    axisBinCenters    = tuple(axes[axisIndex].GetBinCenter(axisBinIndex) for axisIndex, axisBinIndex in enumerate(axisBinIndices))
    yieldsInBin       = {binVarNames[axisIndex] : axisBinCenter for axisIndex, axisBinCenter in enumerate(axisBinCenters)}  # store axis name and bin center
    binIndex          = bins.FindBin(*axisBinCenters)  #!Note! the unpacking works only up to 6 binning dimensions
    fitResultFileName = binFileNames[binIndex].replace("TreeData.root", "ResultsHSMinuit2.root")
    print(f"Reading fit result object 'MinuitResult' from file '{fitResultFileName}'")
    fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
    fitResult     = fitResultFile.Get("MinuitResult")
    fitPars       = fitResult.floatParsFinal()
    for yieldType, yieldParName in YIELD_PAR_NAMES.items():
      yieldParIndex = fitPars.index(yieldParName)
      if yieldParIndex < 0:
        raise IndexError(f"Cannot find parameter '{yieldParName}' in fit parameters {fitPars}")
      yieldPar = fitPars[yieldParIndex]
      yieldsInBin[yieldType] = ufloat(yieldPar.getVal(), yieldPar.getError())  # store yields
    print(f"Read yields for kinematic bin #{binIndex}: {yieldsInBin}")
    yields.append(yieldsInBin)
  return binVarNames, yields


# calculates efficiencies and returns list of dictionaries with efficiencies [ { <binning var> : <bin center>, ..., "Efficiency" : <value> }, ... ]
def calculateEfficiencies(
  binVarNames,  # dict of tuples { <dataset> : (<binning var>, ... ), ... }
  yields        # dict of list of dict { <dataset:> [ { <binning var> : <bin center>, ..., "Signal" : <yield>, "Background" : <yield> }, ... ], ... }
):
  assert binVarNames["Total"] == binVarNames["Found"] == binVarNames["Missing"], "Datasets have different kind of kinematic bins"
  assert len(yields["Total"]) == len(yields["Found"]) == len(yields["Missing"]), "Datasets have different number of kinematic bins"
  efficiencies = []
  for binIndex, yieldFound in enumerate(yields["Found"]):
    yieldMissing = yields["Missing"][binIndex]
    # copy kinematic bin info
    efficiency = {binVarName : yieldFound[binVarName] for binVarName in binVarNames["Found"]}
    # calculate efficiency
    efficiency["Efficiency"] = yieldFound["Signal"] / (yieldFound["Signal"] + yieldMissing["Signal"])
    efficiencies.append(efficiency)
  return efficiencies


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable

  # dataset           = "030730",
  # dataset           = "bggen_2017_01-ver03"
  # dataFileName      = f"../ReactionEfficiency/pippippimpimpmiss_flatTree.{dataset}.root.brufit"
  outputDirName     = "BruFitOutput"
  # kinematicBinnings = [  # list of tuples [ (variable, nmb of bins, min value, max value) ]
  #   ("BeamEnergy", 9, 3.0, 12.0)
  # ]
  dataSets = ["Total", "Found", "Missing"]

  binVarNames = {}
  yields      = {}
  for dataSet in dataSets:
    print(f"Loading yields for '{dataSet}' dataset")
    binVarNames[dataSet], yields[dataSet] = readYieldsFromFitResults(f"{outputDirName}/{dataSet}")
  #TODO implement case when no binning is present

  efficiencies = calculateEfficiencies(binVarNames, yields)
  print(efficiencies)
