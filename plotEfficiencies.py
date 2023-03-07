#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import array
import itertools
import numpy as np

from uncertainties import ufloat

import ROOT

import makePlots  # defines helper functions to generate histograms from data trees

makePlots.setupPlotStyle()


YIELD_PAR_NAMES = {
  "Signal"     : "Yld_SigPdf",
  "Background" : "Yld_BkgPdf"
}

BINNING_VAR_PLOT_INFO = {
  "BeamEnergy"         : {"label" : "E_{beam}",                      "unit" : "GeV"},
  "MissingProtonP"     : {"label" : "#it{p}^{miss}_{kin. fit}",      "unit" : "GeV/c"},
  "MissingProtonTheta" : {"label" : "#it{#theta}^{miss}_{kin. fit}", "unit" : "deg"},
  "MissingProtonPhi"   : {"label" : "#it{#phi}^{miss}_{kin. fit}",   "unit" : "deg"}
}


# reads yields from fit results and returns
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


# calculates efficiencies and returns list of dict with efficiencies [ { <binning var> : <bin center>, ..., "Efficiency" : <value> }, ... ]
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
    print(f"Effiency = {efficiency}")
    efficiencies.append(efficiency)
  return efficiencies


# plot efficiency for 1-dimensional binning
def plotEfficiencies1D(
  binVarNames,   # dict of tuples { <dataset> : (<binning var>, ... ), ... }
  efficiencies,  # list of dict [ { <binning var> : <bin center>, ..., "Efficiency" : <value> }, ... ]
  pdfFileNameSuffix = "",
  particle          = "Proton",
  channel           = "4pi",
  markerSize        = 0.75
):
  assert len(binVarNames["Found"]) == 1, f"This function cannot plot binning with {len(binVarNames['Found'])} variables"
  binVarName  = binVarNames["Found"][0]
  assert binVarName in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binVarName}'"
  binVarLabel = BINNING_VAR_PLOT_INFO[binVarName]["label"]
  binVarUnit  = BINNING_VAR_PLOT_INFO[binVarName]["unit"]
  graphVals = [(efficiency[binVarName], efficiency["Efficiency"]) for efficiency in efficiencies]
  xVals = array.array('d', [graphVal[0]               for graphVal in graphVals])
  yVals = array.array('d', [graphVal[1].nominal_value for graphVal in graphVals])
  yErrs = array.array('d', [graphVal[1].std_dev       for graphVal in graphVals])
  # print(xVals, yVals, yErrs)
  efficienciesKinBinsGraph = ROOT.TGraphErrors(len(graphVals), xVals, yVals, ROOT.nullptr, yErrs)
  efficienciesKinBinsGraph.SetTitle(f"{particle} Track-Finding Efficiency ({channel})")
  efficienciesKinBinsGraph.SetMarkerStyle(ROOT.kFullCircle)
  efficienciesKinBinsGraph.SetMarkerSize(markerSize)
  efficienciesKinBinsGraph.GetXaxis().SetTitle(f"{binVarLabel} ({binVarUnit})")
  efficienciesKinBinsGraph.GetYaxis().SetTitle("Efficiency")
  efficienciesKinBinsGraph.SetMinimum(0)
  efficienciesKinBinsGraph.SetMaximum(1)
  canv = ROOT.TCanvas(f"{particle}_{channel}_mm2_eff_{binVarName}{pdfFileNameSuffix}", "")
  efficienciesKinBinsGraph.Draw("AP")
  #TODO? # indicate value from fit of overall distributions
  line = ROOT.TLine()
  # line.SetLineStyle(ROOT.kDashed)
  # line.DrawLine(efficienciesKinBinsGraph.GetXaxis().GetXmin(), overallEff.nominal_value, efficienciesKinBinsGraph.GetXaxis().GetXmax(), overallEff.nominal_value)
  # indicate weighted average of efficiencies in kinematic bins
  meanEff = np.average(yVals, weights = [1 / (yErr**2) for yErr in yErrs])
  line.SetLineColor(ROOT.kRed + 1)
  line.DrawLine(efficienciesKinBinsGraph.GetXaxis().GetXmin(), meanEff, efficienciesKinBinsGraph.GetXaxis().GetXmax(), meanEff)
  canv.SaveAs(f"{outputDirName}/{canv.GetName()}.pdf")


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable

  outputDirName = "BruFitOutput"
  dataSets      = ["Total", "Found", "Missing"]

  binVarNames = {}
  yields      = {}
  for dataSet in dataSets:
    print(f"Loading yields for '{dataSet}' dataset")
    binVarNames[dataSet], yields[dataSet] = readYieldsFromFitResults(f"{outputDirName}/{dataSet}")
  #TODO implement case when no binning is present

  efficiencies = calculateEfficiencies(binVarNames, yields)

  plotEfficiencies1D(binVarNames, efficiencies)
