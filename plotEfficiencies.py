#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import argparse
import functools
import math
import numpy as np
import os
import sys

from uncertainties import ufloat

import ROOT

import makePlots
import plotFitResults


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


YIELD_PAR_NAMES = {
  "Signal"     : "Yld_SigPdf",
  "Background" : "Yld_BkgPdf"
}


#TODO use dataclasses to store information
#TODO add type annotations
#TODO add docstrings
# reads data histogram in given file and returns its integral
def readDataIntegralFromFitFile(
  fitResultFileName,
  fitVariable,
  binName = ""
):
  print(f"Reading fitted data histogram from file '{fitResultFileName}'")
  fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
  canvName = f"{binName}_{fitVariable}"
  canv = fitResultFile.Get(canvName)
  # get data histogram
  dataFitPad = canv.GetListOfPrimitives().FindObject(f"{canvName}_1")
  dataHist = dataFitPad.GetListOfPrimitives().FindObject("h_DataEvents")
  #TODO reimplement in numpy
  yVals = dataHist.GetY()
  integral = 0
  for i in range(dataHist.GetN()):
    integral += yVals[i]
  fitResultFile.Close()
  return ufloat(integral, math.sqrt(integral))


# reads integrals from data histograms in given directory and returns
#     list of dicts with yields [ { <binning var> : <bin center>, ..., "Signal" : <integral> }, ... ]
#     tuple with binning variables (<binning var>, ... )
def readDataIntegralsFromFitDir(
  fitResultDirName,
  fitVariable
):
  integrals = []
  # read overall yields from fit-result file
  fitResultFileName = f"{fitResultDirName}/ResultsHSMinuit2.root"
  if os.path.isfile(fitResultFileName):
    integralOverall = {"Signal" : readDataIntegralFromFitFile(fitResultFileName, fitVariable)}
    integralOverall.update({"Overall" : None})  # tag overall yields with special axis name
    print(f"Read overall integral: {integralOverall}")
    integrals.append(integralOverall)  # first entry
  # get binning from file
  binningInfo = plotFitResults.getBinningFromDir(fitResultDirName)
  if binningInfo is not None and binningInfo.infos:
    for binInfo in binningInfo.infos:
      integralInBin = {"Signal" : readDataIntegralFromFitFile(binInfo.fitResultFileName, fitVariable, binInfo.name)}
      # copy axis name and bin center
      integralInBin.update(binInfo.centers)
      print(f"Read integral for kinematic bin: {integralInBin}")
      integrals.append(integralInBin)
  return integrals


#TODO reduce code doubling
# reads yields from fit results in given directory and returns
#     list of dicts with yields [ { <binning var> : <bin center>, ..., "Signal" : <yield>, "Background" : <yield> }, ... ]
#     tuple with binning variables (<binning var>, ... )
def readYieldsFromFitDir(fitResultDirName):
  yields = []
  # read overall yields from fit-result file
  fitResultFileName = f"{fitResultDirName}/ResultsHSMinuit2.root"
  if os.path.isfile(fitResultFileName):
    yieldsOverall, _ = plotFitResults.readParValuesFromFitFile(fitResultFileName, YIELD_PAR_NAMES)
    yieldsOverall.update({"Overall" : None})  # tag overall yields with special axis name
    print(f"Read overall yields: {yieldsOverall}")
    yields.append(yieldsOverall)  # first entry
  # get binning from file
  binningInfo = plotFitResults.getBinningFromDir(fitResultDirName)
  if binningInfo is not None and binningInfo.infos:
    for binInfo in binningInfo.infos:
      yieldsInBin, _ = plotFitResults.readParValuesFromFitFile(binInfo.fitResultFileName, YIELD_PAR_NAMES)
      # copy axis name and bin center
      yieldsInBin.update(binInfo.centers)
      print(f"Read yields for kinematic bin: {yieldsInBin}")
      yields.append(yieldsInBin)
  return yields


# calculates efficiencies and returns list of dict with efficiencies [ { <binning var> : <bin center>, ..., "Efficiency" : <value> }, ... ]
def calculateEfficiencies(
  yields,      # dict of lists of dicts with yields { <dataset> : [ { <binning var> : <bin center>, ..., "Signal" : <yield>, "Background" : <yield> }, ... ], ... }
  binVarNames  # list of lists with binning variables for each binning [ [ <binning var>, ... ], ... ]
):
  assert len(yields["Total"]) == len(yields["Found"]) == len(yields["Missing"]), "Datasets have different number of kinematic bins"
  efficiencies = []
  for index, yieldFound in enumerate(yields["Found"]):
    yieldMissing = yields["Missing"][index]
    efficiency = {}
    # calculate efficiency
    efficiency["Efficiency"] = yieldFound["Signal"] / (yieldFound["Signal"] + yieldMissing["Signal"])
    # add axis name and bin center
    efficiency.update({key : yieldFound[key] for key in yieldFound if (key != "Signal") and (key != "Background")})
    print(f"Effiency = {efficiency}")
    efficiencies.append(efficiency)
  return efficiencies


# plots efficiency for 1-dimensional binning
def plotEfficiencies1D(
  efficiencies,  # list of dicts with efficiencies [ { <binning var> : <bin center>, ..., "Efficiency" : <value> }, ... ]
  binningVars,   # tuple with binning variables (<binning var>, ... )
  pdfDirName,    # directory name the PDF file will be written to
  pdfFileNameSuffix = "",
  particle          = "Proton",
  channel           = "4pi",
  markerSize        = 0.75
):
  binningVar = binningVars[0]
  print(f"Plotting efficiency as a function of binning variable '{binningVar}'")
  xVals, yVals, yErrs, binVarLabel, binVarUnit = plotFitResults.getDataPointArrays1D(binningVar, "Efficiency", efficiencies)
  # # set uncertainties to zero as long as they are not estimated well
  # for i in range(len(yErrs)):
  #   yErrs[i] = 0
  # print(xVals, yVals, yErrs)
  efficiencyGraph = ROOT.TGraphErrors(len(xVals), xVals, yVals, ROOT.nullptr, yErrs)
  efficiencyGraph.SetTitle(f"{particle} Track-Finding Efficiency ({channel})")
  efficiencyGraph.SetMarkerStyle(ROOT.kFullCircle)
  efficiencyGraph.SetMarkerSize(markerSize)
  efficiencyGraph.GetXaxis().SetTitle(f"{binVarLabel} ({binVarUnit})")
  efficiencyGraph.GetYaxis().SetTitle("Efficiency")
  efficiencyGraph.SetMinimum(0)
  efficiencyGraph.SetMaximum(1)
  canv = ROOT.TCanvas(f"{particle}_{channel}_mm2_eff_{binningVar}{pdfFileNameSuffix}", "")
  efficiencyGraph.Draw("AP")
  #TODO? # indicate value from fit of overall distributions
  # line = ROOT.TLine()
  # line.SetLineStyle(ROOT.kDashed)
  # line.DrawLine(efficienciesKinBinsGraph.GetXaxis().GetXmin(), overallEff.nominal_value, efficienciesKinBinsGraph.GetXaxis().GetXmax(), overallEff.nominal_value)
  # # indicate weighted average of efficiencies in kinematic bins
  # meanEff = np.average(yVals, weights = [1 / (yErr**2) for yErr in yErrs])
  # line.SetLineColor(ROOT.kRed + 1)
  # line.DrawLine(efficiencyGraph.GetXaxis().GetXmin(), meanEff, efficiencyGraph.GetXaxis().GetXmax(), meanEff)
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  makePlots.setupPlotStyle()
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Plots BruFit results.")
  parser.add_argument("outputDirName", type = str, nargs = "?", default = "./BruFitOutput", help = "The path to the BruFit output directory; (default: '%(default)s')")
  args = parser.parse_args()
  dataSets    = ["Total", "Found", "Missing"]
  fitVariable = "MissingMassSquared_Measured"  #TODO get this info from BruFit's ROOT.Setup

  print("Calculating efficiencies from integrals of data histograms")
  integrals   = {}  # dict of lists of dicts with integrals of data histogram { <dataset:> [ { <binning var> : <bin center>, ..., "Signal" : <integral> }, ... ], ... }
  binVarNames = None  # list of lists with binning variables for each binning [ [ <binning var>, ... ], ... ]
  for dataSet in dataSets:
    print(f"Reading integrals of data histograms for '{dataSet}' dataset")
    fitResultDirName = f"{args.outputDirName}/{dataSet}"
    integrals[dataSet] = readDataIntegralsFromFitDir(fitResultDirName, fitVariable)  # overall integral values
    binVarNamesInDataSet = []
    for binningInfo in plotFitResults.getBinningInfosFromDir(fitResultDirName):
      if binningInfo:
        binVarNamesInDataSet.append(binningInfo.varNames)
        if binningInfo.dirName:
          integrals[dataSet][len(integrals[dataSet]):] = readDataIntegralsFromFitDir(binningInfo.dirName, fitVariable)  # append integral values
    if binVarNames is None:
      binVarNames = binVarNamesInDataSet
    else:
      assert binVarNamesInDataSet == binVarNames, f"The binning variables {binVarNamesInDataSet} of dataset '{dataSet}' are different from the binning variables {binVarNames} of the previous one"
  efficiencies = calculateEfficiencies(integrals, binVarNames)
  if efficiencies:
    if "Overall" in efficiencies[0]:
      print(f"Overall efficiency from data-histogram integral is {efficiencies[0]['Efficiency']}")
    if binVarNames:
      for binningVars in binVarNames:
        if len(binningVars) == 1:
          plotEfficiencies1D(efficiencies, binningVars, args.outputDirName, "_integral")

  #TODO reduce code doubling
  print("Calculating efficiencies from fit results")
  yields      = {}  # dict of lists of dicts with yields { <dataset:> [ { <binning var> : <bin center>, ..., "Signal" : <yield>, "Background" : <yield> }, ... ], ... }
  binVarNames = None  # list of lists with binning variables for each binning [ [ <binning var>, ... ], ... ]
  for dataSet in dataSets:
    print(f"Reading yields for '{dataSet}' dataset")
    fitResultDirName = f"{args.outputDirName}/{dataSet}"
    yields[dataSet] = readYieldsFromFitDir(fitResultDirName)  # overall yield values
    binVarNamesInDataSet = []
    for binningInfo in plotFitResults.getBinningInfosFromDir(fitResultDirName):
      if binningInfo:
        binVarNamesInDataSet.append(binningInfo.varNames)
        if binningInfo.dirName:
          yields[dataSet][len(yields[dataSet]):] = readYieldsFromFitDir(binningInfo.dirName)  # append yield values
    if binVarNames is None:
      binVarNames = binVarNamesInDataSet
    else:
      assert binVarNamesInDataSet == binVarNames, f"The binning variables {binVarNamesInDataSet} of dataset '{dataSet}' are different from the binning variables {binVarNames} of the previous one"
  efficiencies = calculateEfficiencies(yields, binVarNames)
  if efficiencies:
    if "Overall" in efficiencies[0]:
      print(f"Overall efficiency from fits is {efficiencies[0]['Efficiency']}")
    if binVarNames:
      for binningVars in binVarNames:
        if len(binningVars) == 1:
          plotEfficiencies1D(efficiencies, binningVars, args.outputDirName)
