#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import array
import os

import ROOT
ROOT.PyConfig.DisableRootLogon = True

from uncertainties import ufloat

import plotEfficiencies


DATASET_COLORS = {
  "Total"   : ROOT.kBlack,
  "Found"   : ROOT.kGreen + 2,
  "Missing" : ROOT.kRed + 1
}


# plots fit results for kinematic bins
def plotFitResults(bins):
  # assume that lists returned by ROOT.Bins.GetBinNames() and ROOT.Bins.GetFileNames() have the same ordering
  binNames = [str(binName) for binName in bins.GetBinNames()]
  fitResultFileNames = [str(fileName).replace("TreeData.root", "ResultsHSMinuit2.root") for fileName in bins.GetFileNames()]
  for index, fitResultFileName  in enumerate(fitResultFileNames):
    fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
    print(f"Plotting fit result in '{fitResultFileName}'")
    canvName = f"{binNames[index]}_{fitVariable}"
    canv = fitResultFile.Get(canvName)
    # remove TPaveText with fit parameters
    dataFitPad = canv.GetListOfPrimitives().FindObject(f"{canvName}_1")
    paramBox = dataFitPad.GetListOfPrimitives().FindObject(f"{binNames[index]}TotalPDF_paramBox")
    dataFitPad.GetListOfPrimitives().Remove(paramBox)
    # paramBox.SetBorderSize(0)
    # paramBox.SetFillStyle(0)
    canv.SaveAs(f"{os.path.dirname(fitResultFileName)}/{canv.GetName()}.pdf")
    fitResultFile.Close()


# reads fit result in given file and returns
#    dict with parameter values { <par name> : <par value>, ... }
#    tuple with parameter names
def readParValuesFromFitFile(fitResultFileName):
  print(f"Reading fit result object 'MinuitResult' from file '{fitResultFileName}'")
  fitResultFile  = ROOT.TFile.Open(fitResultFileName, "READ")
  fitResult      = fitResultFile.Get("MinuitResult")
  fitPars        = fitResult.floatParsFinal()
  parValuesInBin = {}
  for fitParIndex in range(fitPars.getSize()):
    fitPar = fitPars[fitParIndex]
    parValuesInBin[fitPar.GetName()] = ufloat(fitPar.getVal(), fitPar.getError())
  fitResultFile.Close()
  parNamesInBin = tuple(parValuesInBin.keys())
  return parValuesInBin, parNamesInBin


#!TODO unify this and the above function with readYieldsFromFitDir() in plotEfficiencies
# reads parameter values from fit results in given directory and returns
#     tuple with binning variables (<binning var>, ... )
#     list of dicts with parameter values [ { <binning var> : <bin center>, ..., <par name> : <par value>, ... }, ... ]
def readParValuesFromFitDir(bins, binVarNames):
  parValues = []
  parNames = None  # used to compare parameter names in current and previous bin
  for fitResultFileName in plotEfficiencies.getFitResultFileNames(bins, binVarNames):
    parValuesInBin, parNamesInBin = readParValuesFromFitFile(fitResultFileName["FitResultFileName"])
    # ensure that the parameters set of the fit function is the same in all kinematic bins
    if parNames is not None:
      assert parNamesInBin == parNames, f"The parameter set {parNamesInBin} of this bin is different from the parameter set {parNames} of the previous one"
    parNames = parNamesInBin
    # copy axis name and bin center
    parValuesInBin.update({key : fitResultFileName[key] for key in fitResultFileName if key != "FitResultFileName"})
    print(f"Read parameter values for kinematic bin: {parValuesInBin}")
    parValues.append(parValuesInBin)
  return parValues, parNames


# helper function that draws zero line when needed
def drawZeroLine(xAxis, yAxis, style = ROOT.kDashed, color = ROOT.kBlack):
  if (yAxis.GetXmin() < 0) and (yAxis.GetXmax() > 0):
    zeroLine = ROOT.TLine()
    zeroLine.SetLineStyle(style)
    zeroLine.SetLineColor(color)
    return zeroLine.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 0, xAxis.GetBinUpEdge(xAxis.GetLast()), 0)


# overlay parameter value for datasets for 1-dimensional binning
def plotParValue1D(
  parValues,      # dict of lists of dicts with parameter values { <dataset> : [ { <binning var> : <bin center>, ..., <par name> : <par value>, ... }, ... ], ... }
  parName,        # name of parameter to plot
  binVarNames,    # tuple with binning variables (<binning var>, ... )
  outputDirName,  # directory name the PDF file will be written to
  pdfFileNameSuffix = "",
  particle          = "Proton",
  channel           = "4pi",
  markerSize        = 0.75
):
  binVarName  = binVarNames[0]
  print(f"Plotting parameter '{parName}' as a function of binning variable '{binVarName}'")
  parValueMultiGraph = ROOT.TMultiGraph()
  parValueGraphs = {}  # store graphs here to keep them in memory
  for dataSet in parValues:
    xVals, yVals, yErrs, binVarLabel, binVarUnit = plotEfficiencies.getDataPointArrays1D(binVarName, parName, parValues[dataSet])
    # print(dataSet, xVals, yVals, yErrs, binVarLabel, binVarUnit)
    graph = parValueGraphs[dataSet] = ROOT.TGraphErrors(len(xVals), xVals, yVals, ROOT.nullptr, yErrs)
    graph.SetTitle(dataSet)
    graph.SetMarkerStyle(ROOT.kCircle)
    graph.SetMarkerSize(markerSize)
    graph.SetMarkerColor(DATASET_COLORS[dataSet])
    graph.SetLineColor(DATASET_COLORS[dataSet])
    parValueMultiGraph.Add(graph)
  parValueMultiGraph.SetTitle(f"Fit parameter {parName}, {particle} ({channel})")
  parValueMultiGraph.GetXaxis().SetTitle(f"{binVarLabel} ({binVarUnit})")
  parValueMultiGraph.GetYaxis().SetTitle(parName)
  canv = ROOT.TCanvas(f"{particle}_{channel}_{parName}_{binVarName}{pdfFileNameSuffix}", "")
  parValueMultiGraph.Draw("APZ")
  canv.BuildLegend()
  drawZeroLine(parValueMultiGraph.GetXaxis(), parValueMultiGraph.GetYaxis())
  canv.SaveAs(f"{outputDirName}/{canv.GetName()}.pdf")


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable

  outputDirName = "BruFitOutput"
  dataSets      = ["Total", "Found", "Missing"]
  fitVariable   = "MissingMassSquared_Measured"  #TODO read this from ROOT.Setup

  parValues   = {}
  parNames    = None
  binVarNames = None
  for dataSet in dataSets:
    bins, binVarNamesInDataSet = plotEfficiencies.getBinningFromFile(f"{outputDirName}/{dataSet}")
    if binVarNames is not None:
      assert binVarNamesInDataSet == binVarNames, f"The binning variables {binVarNamesInDataSet} of dataset '{dataSet}' are different from the binning variables {binVarNames} of the previous one"
    binVarNames = binVarNamesInDataSet
    if bins is not None:
      plotFitResults(bins)
      parValues[dataSet], parNamesInDataSet = readParValuesFromFitDir(bins, binVarNames)
      if parNames is not None:
        assert parNamesInDataSet == parNames, f"The parameter set {parNamesInDataSet} of dataset '{dataSet}' is different from the parameter set {parNames} of the previous one"
      parNames = parNamesInDataSet
  for parName in parNames:
    plotParValue1D(parValues, parName, binVarNames, outputDirName)
