#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import array
import itertools
import numpy as np
import os

from uncertainties import ufloat

import ROOT

import makePlots


YIELD_PAR_NAMES = {
  "Signal"     : "Yld_SigPdf",
  "Background" : "Yld_BkgPdf"
}

#TODO move to plotFitResults
BINNING_VAR_PLOT_INFO = {
  "BeamEnergy"         : {"label" : "E_{beam}",                      "unit" : "GeV"},
  "MissingProtonP"     : {"label" : "#it{p}^{miss}_{kin. fit}",      "unit" : "GeV/c"},
  "MissingProtonTheta" : {"label" : "#it{#theta}^{miss}_{kin. fit}", "unit" : "deg"},
  "MissingProtonPhi"   : {"label" : "#it{#phi}^{miss}_{kin. fit}",   "unit" : "deg"}
}


#TODO move to plotFitResults
# returns
#     ROOT.Bins object in from binning file in given directory
#     tuple with binning variables (<binning var>, ... )
def getBinningFromFile(fitResultDirName):
  binningFileName = f"{fitResultDirName}/DataBinsConfig.root"
  if not os.path.isfile(binningFileName):
    return None, None
  print(f"Loading binning from file '{binningFileName}'")
  bins = ROOT.Bins("HSBins", binningFileName)
  print("Found binning:")
  bins.PrintAxis()
  axes        = bins.GetVarAxis()
  binVarNames = tuple(axis.GetName() for axis in axes)
  return bins, binVarNames


#TODO move to plotFitResults
# returns list of dicts with file names of fit results [ { <binning var> : <bin center>, ..., "FitResultFileName" : <name> }, ... ]
def getFitResultFileNames(
  bins,
  binVarNames  # tuple with binning variables (<binning var>, ... )
):
  fitResultFileNames = []
  axes         = bins.GetVarAxis()
  binFileNames = [str(fileName) for fileName in bins.GetFileNames()]
  axisBinIndexRanges = tuple(range(1, axis.GetNbins() + 1) for axis in axes)
  for axisBinIndices in itertools.product(*axisBinIndexRanges):  # loop over all tuples of bin indices for the axes
    axisBinCenters         = tuple(axes[axisIndex].GetBinCenter(axisBinIndex) for axisIndex, axisBinIndex in enumerate(axisBinIndices))
    binIndex               = bins.FindBin(*axisBinCenters)  #!Note! the unpacking works only for up to 6 binning dimensions
    fitResultFileNameInBin = {"FitResultFileName" : binFileNames[binIndex].replace("TreeData.root", "ResultsHSMinuit2.root")}
    fitResultFileNameInBin.update({binVarNames[axisIndex] : axisBinCenter for axisIndex, axisBinCenter in enumerate(axisBinCenters)})
    fitResultFileNames.append(fitResultFileNameInBin)
  return fitResultFileNames


#TODO move to plotFitResults
# reads yields from fit result in given file and returns dict with yields { "Signal" : <yield>, "Background" : <yield> }
def readYieldsFromFitFile(fitResultFileName):
  print(f"Reading fit result object 'MinuitResult' from file '{fitResultFileName}'")
  fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
  fitResult     = fitResultFile.Get("MinuitResult")
  fitPars       = fitResult.floatParsFinal()
  yieldsInBin   = {}
  for yieldType, yieldParName in YIELD_PAR_NAMES.items():
    yieldParIndex = fitPars.index(yieldParName)
    if yieldParIndex < 0:
      raise IndexError(f"Cannot find parameter '{yieldParName}' in fit parameters {fitPars} in file '{fitResultFileName}'")
    yieldPar = fitPars[yieldParIndex]
    yieldsInBin[yieldType] = ufloat(yieldPar.getVal(), yieldPar.getError())  # store yields
  fitResultFile.Close()
  return yieldsInBin


#TODO move to plotFitResults
# reads yields from fit results in given directory and returns
#     list of dicts with yields [ { <binning var> : <bin center>, ..., "Signal" : <yield>, "Background" : <yield> }, ... ]
#     tuple with binning variables (<binning var>, ... )
def readYieldsFromFitDir(fitResultDirName):
  yields = []
  # read overall yields from fit-result file
  fitResultFileName = f"{fitResultDirName}/ResultsHSMinuit2.root"
  if os.path.isfile(fitResultFileName):
    yieldsInBin = readYieldsFromFitFile(fitResultFileName)
    yieldsInBin.update({"Overall" : None})  # tag overall yields with special axis name
    print(f"Read overall yields: {yieldsInBin}")
    yields.append(yieldsInBin)  # first entry
  # get binning from file
  bins, binVarNames = getBinningFromFile(fitResultDirName)
  if bins is not None:
    for fitResultFileName in getFitResultFileNames(bins, binVarNames):
      yieldsInBin = readYieldsFromFitFile(fitResultFileName["FitResultFileName"])
      # copy axis name and bin center
      yieldsInBin.update({key : fitResultFileName[key] for key in fitResultFileName if key != "FitResultFileName"})
      print(f"Read yields for kinematic bin: {yieldsInBin}")
      yields.append(yieldsInBin)
  return yields, binVarNames


# calculates efficiencies and returns list of dict with efficiencies [ { <binning var> : <bin center>, ..., "Efficiency" : <value> }, ... ]
def calculateEfficiencies(
  yields,      # dict of lists of dicts with yields { <dataset:> [ { <binning var> : <bin center>, ..., "Signal" : <yield>, "Background" : <yield> }, ... ], ... }
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


#TODO move to plotFitResults
# returns
#     arrays of x, y, and y-uncertainty values
#     x-axis label and unit
def getDataPointArrays1D(
  binVarName,  # name of the binning variable, i.e. x-axis
  valueName,   # name of value, i.e. y-axis
  values       # list of dicts with data values [ { <binning var> : <bin center>, ..., <value name> : <value> }, ... ]
):
  assert binVarName in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binVarName}'"
  binVarLabel = BINNING_VAR_PLOT_INFO[binVarName]["label"]
  binVarUnit  = BINNING_VAR_PLOT_INFO[binVarName]["unit"]
  graphVals = [(value[binVarName], value[valueName]) for value in values if (binVarName in value) and (valueName in value)]
  if not graphVals:
    print("No data to plot")
    return
  xVals = array.array('d', [graphVal[0]               for graphVal in graphVals])
  yVals = array.array('d', [graphVal[1].nominal_value for graphVal in graphVals])
  yErrs = array.array('d', [graphVal[1].std_dev       for graphVal in graphVals])
  return xVals, yVals, yErrs, binVarLabel, binVarUnit


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
  xVals, yVals, yErrs, binVarLabel, binVarUnit = getDataPointArrays1D(binVarName, "Efficiency", efficiencies)
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
  # indicate weighted average of efficiencies in kinematic bins
  meanEff = np.average(yVals, weights = [1 / (yErr**2) for yErr in yErrs])
  line.SetLineColor(ROOT.kRed + 1)
  line.DrawLine(efficiencyGraph.GetXaxis().GetXmin(), meanEff, efficiencyGraph.GetXaxis().GetXmax(), meanEff)
  canv.SaveAs(f"{outputDirName}/{canv.GetName()}.pdf")


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable
  makePlots.setupPlotStyle()

  outputDirName = "BruFitOutput"
  dataSets      = ["Total", "Found", "Missing"]

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
