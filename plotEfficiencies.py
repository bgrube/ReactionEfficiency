#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import math
import numpy as np
import os

from uncertainties import ufloat

import ROOT

import makePlots
import plotFitResults


YIELD_PAR_NAMES = {
  "Signal"     : "Yld_SigPdf",
  "Background" : "Yld_BkgPdf"
}


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
  bins, binVarNames = plotFitResults.getBinningFromDir(fitResultDirName)
  if bins is not None:
    for fitResultFileName in plotFitResults.getFitResultFileNames(bins, binVarNames):
      integralInBin = {"Signal" : readDataIntegralFromFitFile(fitResultFileName["FitResultFileName"], fitVariable, fitResultFileName["BinName"])}
      # copy axis name and bin center
      integralInBin.update({key : fitResultFileName[key] for key in fitResultFileName if key != "FitResultFileName"})
      print(f"Read integral for kinematic bin: {integralInBin}")
      integrals.append(integralInBin)
  return integrals, binVarNames


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
  bins, binVarNames = plotFitResults.getBinningFromDir(fitResultDirName)
  if bins is not None:
    for fitResultFileName in plotFitResults.getFitResultFileNames(bins, binVarNames):
      yieldsInBin, _ = plotFitResults.readParValuesFromFitFile(fitResultFileName["FitResultFileName"], YIELD_PAR_NAMES)
      # copy axis name and bin center
      yieldsInBin.update({key : fitResultFileName[key] for key in fitResultFileName if key != "FitResultFileName"})
      print(f"Read yields for kinematic bin: {yieldsInBin}")
      yields.append(yieldsInBin)
  return yields, binVarNames


# calculates efficiencies and returns list of dict with efficiencies [ { <binning var> : <bin center>, ..., "Efficiency" : <value> }, ... ]
def calculateEfficiencies(
  yields,      # dict of lists of dicts with yields { <dataset> : [ { <binning var> : <bin center>, ..., "Signal" : <yield>, "Background" : <yield> }, ... ], ... }
  binVarNames  # dict of tuples with binning variables { <dataset> : (<binning var>, ... ), ... }
):
  assert binVarNames["Total"] == binVarNames["Found"] == binVarNames["Missing"], "Datasets have different kind of kinematic bins"
  assert len(yields["Total"]) == len(yields["Found"]) == len(yields["Missing"]), "Datasets have different number of kinematic bins"
  efficiencies = []
  for index, yieldFound in enumerate(yields["Found"]):
    yieldMissing = yields["Missing"][index]
    efficiency = {}
    # calculate efficiency
    efficiency["Efficiency"] = yieldFound["Signal"] / (yieldFound["Signal"] + yieldMissing["Signal"])
    # add axis name and bin center
    if (index == 0) and ("Overall" in yieldFound):
      efficiency.update({"Overall" : None})  # tag overall efficiency with special axis name
    elif binVarNames:
      efficiency.update({binVarName : yieldFound[binVarName] for binVarName in binVarNames["Found"]})
    print(f"Effiency = {efficiency}")
    efficiencies.append(efficiency)
  return efficiencies


# plots efficiency for 1-dimensional binning
def plotEfficiencies1D(
  efficiencies,   # list of dicts with efficiencies [ { <binning var> : <bin center>, ..., "Efficiency" : <value> }, ... ]
  binVarNames,    # dict of tuples with binning variables { <dataset> : (<binning var>, ... ), ... }
  outputDirName,  # directory name the PDF file will be written to
  pdfFileNameSuffix = "",
  particle          = "Proton",
  channel           = "4pi",
  markerSize        = 0.75
):
  binVarName  = binVarNames["Found"][0]
  print(f"Plotting efficiency as a function of binning variable '{binVarName}'")
  xVals, yVals, yErrs, binVarLabel, binVarUnit = plotFitResults.getDataPointArrays1D(binVarName, "Efficiency", efficiencies)
  # set uncertainties to zero as long as they are not estimated well
  for i in range(len(yErrs)):
    yErrs[i] = 0
  # print(xVals, yVals, yErrs)
  efficiencyGraph = ROOT.TGraphErrors(len(xVals), xVals, yVals, ROOT.nullptr, yErrs)
  efficiencyGraph.SetTitle(f"{particle} Track-Finding Efficiency ({channel})")
  efficiencyGraph.SetMarkerStyle(ROOT.kFullCircle)
  efficiencyGraph.SetMarkerSize(markerSize)
  efficiencyGraph.GetXaxis().SetTitle(f"{binVarLabel} ({binVarUnit})")
  efficiencyGraph.GetYaxis().SetTitle("Efficiency")
  efficiencyGraph.SetMinimum(0)
  efficiencyGraph.SetMaximum(1)
  canv = ROOT.TCanvas(f"{particle}_{channel}_mm2_eff_{binVarName}{pdfFileNameSuffix}", "")
  efficiencyGraph.Draw("AP")
  #TODO? # indicate value from fit of overall distributions
  line = ROOT.TLine()
  # line.SetLineStyle(ROOT.kDashed)
  # line.DrawLine(efficienciesKinBinsGraph.GetXaxis().GetXmin(), overallEff.nominal_value, efficienciesKinBinsGraph.GetXaxis().GetXmax(), overallEff.nominal_value)
  # # indicate weighted average of efficiencies in kinematic bins
  # meanEff = np.average(yVals, weights = [1 / (yErr**2) for yErr in yErrs])
  # line.SetLineColor(ROOT.kRed + 1)
  # line.DrawLine(efficiencyGraph.GetXaxis().GetXmin(), meanEff, efficiencyGraph.GetXaxis().GetXmax(), meanEff)
  canv.SaveAs(f"{outputDirName}/{canv.GetName()}.pdf")


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  makePlots.setupPlotStyle()
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable

  outputDirName = "./BruFitOutput"
  dataSets      = ["Total", "Found", "Missing"]
  fitVariable   = "MissingMassSquared_Measured"  #TODO get this info from BruFit's ROOT.Setup

  print("Calculating efficiencies from integrals of data histograms")
  integrals   = {}  # dict of lists of dicts with integrals of data histogram { <dataset:> [ { <binning var> : <bin center>, ..., "Signal" : <integral> }, ... ], ... }
  binVarNames = {}  # dict of tuples with binning variables { <dataset> : (<binning var>, ... ), ... }
  for dataSet in dataSets:
    print(f"Reading integrals of data histograms for '{dataSet}' dataset")
    integrals[dataSet], binVarNames[dataSet] = readDataIntegralsFromFitDir(f"{outputDirName}/{dataSet}", fitVariable)

  efficiencies = calculateEfficiencies(integrals, binVarNames)
  if efficiencies:
    if "Overall" in efficiencies[0]:
      print(f"Overall efficiency from data-histogram integral in directory '{outputDirName}' is {efficiencies[0]['Efficiency']}")
    if binVarNames["Found"] and (len(binVarNames["Found"]) == 1):
      plotEfficiencies1D(efficiencies, binVarNames, outputDirName, "_integral")

  print("Calculating efficiencies from fit results")
  yields      = {}  # dict of lists of dicts with yields { <dataset:> [ { <binning var> : <bin center>, ..., "Signal" : <yield>, "Background" : <yield> }, ... ], ... }
  binVarNames = {}  # dict of tuples with binning variables { <dataset> : (<binning var>, ... ), ... }
  for dataSet in dataSets:
    print(f"Loading yields for '{dataSet}' dataset")
    yields[dataSet], binVarNames[dataSet] = readYieldsFromFitDir(f"{outputDirName}/{dataSet}")

  efficiencies = calculateEfficiencies(yields, binVarNames)
  if efficiencies:
    if "Overall" in efficiencies[0]:
      print(f"Overall efficiency for fits in directory '{outputDirName}' is {efficiencies[0]['Efficiency']}")
    if binVarNames["Found"] and (len(binVarNames["Found"]) == 1):
      plotEfficiencies1D(efficiencies, binVarNames, outputDirName)
