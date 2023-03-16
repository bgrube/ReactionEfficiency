#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import numpy as np
import os

import ROOT

import makePlots
import plotFitResults


YIELD_PAR_NAMES = {
  "Signal"     : "Yld_SigPdf",
  "Background" : "Yld_BkgPdf"
}


# reads yields from fit results in given directory and returns
#     list of dicts with yields [ { <binning var> : <bin center>, ..., "Signal" : <yield>, "Background" : <yield> }, ... ]
#     tuple with binning variables (<binning var>, ... )
def readYieldsFromFitDir(fitResultDirName):
  yields = []
  # read overall yields from fit-result file
  fitResultFileName = f"{fitResultDirName}/ResultsHSMinuit2.root"
  if os.path.isfile(fitResultFileName):
    yieldsInBin, _ = plotFitResults.readParValuesFromFitFile(fitResultFileName, YIELD_PAR_NAMES)
    yieldsInBin.update({"Overall" : None})  # tag overall yields with special axis name
    print(f"Read overall yields: {yieldsInBin}")
    yields.append(yieldsInBin)  # first entry
  # get binning from file
  bins, binVarNames = plotFitResults.getBinningFromFile(fitResultDirName)
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

  outputDirName = "./BruFitOutput"
  dataSets      = ["Total", "Found", "Missing"]

  yields      = {}  # dict of lists of dicts with yields { <dataset:> [ { <binning var> : <bin center>, ..., "Signal" : <yield>, "Background" : <yield> }, ... ], ... }
  binVarNames = {}  # dict of tuples with binning variables { <dataset> : (<binning var>, ... ), ... }
  for dataSet in dataSets:
    print(f"Loading yields for '{dataSet}' dataset")
    yields[dataSet], binVarNames[dataSet] = readYieldsFromFitDir(f"{outputDirName}/{dataSet}")

  #TODO calculate efficiencies also from histogram integrals
  efficiencies = calculateEfficiencies(yields, binVarNames)
  if efficiencies:
    if "Overall" in efficiencies[0]:
      print(f"Overall efficiency for fits in directory '{outputDirName}' is {efficiencies[0]['Efficiency']}")
    if binVarNames["Found"] and (len(binVarNames["Found"]) == 1):
      plotEfficiencies1D(efficiencies, binVarNames, outputDirName)
