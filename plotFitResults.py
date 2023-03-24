#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import array
import collections
import ctypes
import itertools
import os

from uncertainties import ufloat

import ROOT
if __name__ == "__main__":
  ROOT.PyConfig.DisableRootLogon = True  # do not change style of canvases loaded from fit result files

import makePlots


BINNING_VAR_PLOT_INFO = {
  "BeamEnergy"         : {"label" : "E_{beam}",                      "unit" : "GeV"},
  "MissingProtonP"     : {"label" : "#it{p}^{miss}_{kin. fit}",      "unit" : "GeV/c"},
  "MissingProtonTheta" : {"label" : "#it{#theta}^{miss}_{kin. fit}", "unit" : "deg"},
  "MissingProtonPhi"   : {"label" : "#it{#phi}^{miss}_{kin. fit}",   "unit" : "deg"}
}

DATASET_COLORS = {
  "Total"   : ROOT.kBlack,
  "Found"   : ROOT.kGreen + 2,
  "Missing" : ROOT.kRed + 1
}

REMOVE_PARAM_BOX = False


# returns
#     ROOT.Bins object in from binning file in given directory
#     tuple with binning variables (<binning var>, ... )
def getBinningFromDir(fitResultDirName):
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


# returns list of tuples with binning info [ ( <ROOT.Bins object>, binning tuple (<binning var>, ... ) ), ... ]
def getBinningsFromDir(fitResultDirName):
  binnings = []
  # find all subdirectories with binning files
  subDirNames = sorted([entry.path for entry in os.scandir(fitResultDirName) if entry.is_dir() and os.path.isfile(f"{entry.path}/DataBinsConfig.root")])
  # get binning info from all subdirectories
  for subDirName in subDirNames:
    print(f"Found binning info in directory '{subDirName}'")
    binnings.append(getBinningFromDir(subDirName))
  return binnings


# returns list of dicts with file names of fit results [ { <binning var> : <bin center>, ..., "FitResultFileName" : <name>, "BinName" : <name> }, ... ]
def getFitResultFileNames(
  bins,
  binVarNames  # tuple with binning variables (<binning var>, ... )
):
  fitResultFileNames = []
  axes         = bins.GetVarAxis()
  # assume that lists returned by ROOT.Bins.GetBinNames() and ROOT.Bins.GetFileNames() have the same ordering
  binNames     = [str(binName) for binName in bins.GetBinNames()]
  binFileNames = [f"{os.path.dirname(str(fileName))}/ResultsHSMinuit2.root" for fileName in bins.GetFileNames()]
  axisBinIndexRanges = tuple(range(1, axis.GetNbins() + 1) for axis in axes)
  for axisBinIndices in itertools.product(*axisBinIndexRanges):  # loop over all tuples of bin indices for the axes
    axisBinCenters         = tuple(axes[axisIndex].GetBinCenter(axisBinIndex) for axisIndex, axisBinIndex in enumerate(axisBinIndices))
    binIndex               = bins.FindBin(*axisBinCenters)  #!Note! the unpacking works only for up to 6 binning dimensions
    fitResultFileNameInBin = {"FitResultFileName" : binFileNames[binIndex], "BinName" : binNames[binIndex]}
    fitResultFileNameInBin.update({binVarNames[axisIndex] : axisBinCenter for axisIndex, axisBinCenter in enumerate(axisBinCenters)})
    fitResultFileNames.append(fitResultFileNameInBin)
  return fitResultFileNames


# reads fit result in given file and returns
#    dict with parameter values { <par name> : <par value>, ... }
#    tuple with parameter names (<par name>, ... )
def readParValuesFromFitFile(
  fitResultFileName,
  fitParNamesToRead = None  # if set to dict { <new par name> : <par name>, ... } only the given parameters are read, where `new par name` is the key used in the output
):
  print(f"Reading fit result object 'MinuitResult' from file '{fitResultFileName}'")
  fitResultFile  = ROOT.TFile.Open(fitResultFileName, "READ")
  fitResult      = fitResultFile.Get("MinuitResult")
  fitPars        = fitResult.floatParsFinal()
  parValuesInBin = {}
  if isinstance(fitParNamesToRead, collections.Mapping):
    # fitParNamesToRead is some kind of dict
    # read only selected fit parameters; and use new key names for them
    for fitParKey, fitParName in fitParNamesToRead.items():
      fitParIndex = fitPars.index(fitParName)
      if fitParIndex < 0:
        #TODO try using printStream() with ROOT.cout to improve output
        print(f"Cannot find parameter '{fitParName}' in fit parameters {fitPars} in file '{fitResultFileName}'. Skipping parameter.")
        continue
      fitPar = fitPars[fitParIndex]
      parValuesInBin[fitParKey] = ufloat(fitPar.getVal(), fitPar.getError())
  else:
    # read all parameters
    parValuesInBin.update({fitPar.GetName() : ufloat(fitPar.getVal(), fitPar.getError()) for fitPar in fitPars})
  parNamesInBin = tuple(parValuesInBin.keys())
  for i in range(fitResult.numStatusHistory()):
    print(f"    Fit result status {i}: {fitResult.statusCodeHistory(i)} = {fitResult.statusLabelHistory(i)}")
  fitResultFile.Close()
  return parValuesInBin, parNamesInBin


# reads parameter values from fit results in given directory and returns
#     list of dicts with parameter values [ { <binning var> : <bin center>, ..., <par name> : <par value>, ... }, ... ]
#     tuple with parameter names (<par name>, ... )
def readParValuesFromFitDir(bins, binVarNames):
  parValues = []
  parNames = None  # used to compare parameter names in current and previous bin
  for fitResultFileName in getFitResultFileNames(bins, binVarNames):
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
def drawZeroLine(obj, style = ROOT.kDashed, color = ROOT.kBlack):
  objType = obj.IsA().GetName()
  if (objType == "TCanvas") or (objType == "TPad"):
    xMin = ctypes.c_double()
    xMax = ctypes.c_double()
    yMin = ctypes.c_double()
    yMax = ctypes.c_double()
    obj.GetRangeAxis(xMin, yMin, xMax, yMax)
    if (yMin.value < 0) and (yMax.value > 0):
      zeroLine = ROOT.TLine()
      zeroLine.SetLineStyle(style)
      zeroLine.SetLineColor(color)
      return zeroLine.DrawLine(xMin, 0, xMax, 0)
  elif objType.startswith("TH") or objType.startswith("TGraph") or objType.startswith("TMulti"):
    xAxis = obj.GetXaxis()
    yAxis = obj.GetYaxis()
    if (yAxis.GetXmin() < 0) and (yAxis.GetXmax() > 0):
      zeroLine = ROOT.TLine()
      zeroLine.SetLineStyle(style)
      zeroLine.SetLineColor(color)
      return zeroLine.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 0, xAxis.GetBinUpEdge(xAxis.GetLast()), 0)
    elif (yAxis.GetXmin() > 0) and (yAxis.GetXmax() > 0):
      obj.SetMinimum(0)
  else:
    raise TypeError(f"drawZeroLine() not (yet) implemented for object of type '{objType}'")


# plots fit result given directory and saves PDF
def plotFitResult(
  fitResultDirName,
  fitVariable,
  binName = "",
  outputDirName = None  # overrides default PDF output path (i.e. same dir as fit result file) if set
):
  fitResultFileName = f"{fitResultDirName}/ResultsHSMinuit2.root"
  if not os.path.isfile(fitResultFileName):
    print(f"Cannot find file '{fitResultFileName}'; skipping")
    return
  print(f"Plotting fit result in file '{fitResultFileName}'")
  fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
  canvName = f"{binName}_{fitVariable}"
  canv = fitResultFile.Get(canvName)
  # improve TPaveText with fit parameters
  dataFitPad = canv.GetListOfPrimitives().FindObject(f"{canvName}_1")
  paramBox = dataFitPad.GetListOfPrimitives().FindObject(f"{binName}TotalPDF_paramBox")
  if REMOVE_PARAM_BOX:
    # remove box completely
    dataFitPad.GetListOfPrimitives().Remove(paramBox)
  else:
    # only remove filled frame
    paramBox.SetBorderSize(0)
    paramBox.SetFillStyle(0)
  drawZeroLine(dataFitPad)
  pdfFileName = ("Overall" if binName == "" else "") + f"{canv.GetName()}.pdf"
  if outputDirName:
    canv.SaveAs(f"{outputDirName}/{pdfFileName}")
  else:
    canv.SaveAs(f"{os.path.dirname(fitResultFileName)}/{pdfFileName}")
  fitResultFile.Close()


# plots fit results for all kinematic bins
def plotFitResults(
  bins,
  fitVariable,
  outputDirName = None
):
  # assume that lists returned by ROOT.Bins.GetBinNames() and ROOT.Bins.GetFileNames() have the same ordering
  binNames = [str(binName) for binName in bins.GetBinNames()]
  fitResultDirNames = [os.path.dirname(str(fileName)) for fileName in bins.GetFileNames()]
  for index, fitResultDirName  in enumerate(fitResultDirNames):
    plotFitResult(fitResultDirName, fitVariable, binNames[index], outputDirName)


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


# overlay parameter value for datasets for 1-dimensional binning
def plotParValue1D(
  parValues,      # dict of lists of dicts with parameter values { <dataset> : [ { <binning var> : <bin center>, ..., <par name> : <par value>, ... }, ... ], ... }
  parName,        # name of parameter to plot
  binningVars,    # tuple with binning variables (<binning var>, ... )
  outputDirName,  # directory name the PDF file will be written to
  pdfFileNameSuffix = "",
  particle          = "Proton",
  channel           = "4pi",
  markerSize        = 0.75
):
  binVarName  = binningVars[0]
  print(f"Plotting parameter '{parName}' as a function of binning variable '{binVarName}'")
  parValueMultiGraph = ROOT.TMultiGraph()
  parValueGraphs = {}  # store graphs here to keep them in memory
  for dataSet in parValues:
    xVals, yVals, yErrs, binVarLabel, binVarUnit = getDataPointArrays1D(binVarName, parName, parValues[dataSet])
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
  drawZeroLine(parValueMultiGraph)
  canv.SaveAs(f"{outputDirName}/{canv.GetName()}.pdf")


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable

  outputDirName = "./BruFitOutput"
  dataSets      = ["Total", "Found", "Missing"]
  fitVariable   = "MissingMassSquared_Measured"  #TODO get this info from BruFit's ROOT.Setup

  parValues   = {}  # dict of lists of dicts with parameter values { <dataset> : [ { <binning var> : <bin center>, ..., <par name> : <par value>, ... }, ... ], ... }
  parNames    = None
  binVarNames = None  # list of lists with binning variables for each binning [ [ <binning var>, ... ], ... ]
  for dataSet in dataSets:
    fitResultDirName  = f"{outputDirName}/{dataSet}"
    print(f"Plotting overall fit result")
    plotFitResult(fitResultDirName, fitVariable)
    binVarNamesInDataSet = []
    for bins, binningVars in getBinningsFromDir(fitResultDirName):
      binVarNamesInDataSet.append(binningVars)
      if bins is not None:
        print(f"Plotting fit results for binning variable(s) '{binningVars}'")
        plotFitResults(bins, fitVariable, fitResultDirName)
        if not dataSet in parValues:
          parValues[dataSet] = []
        parValues[dataSet][len(parValues[dataSet]):], parNamesInBinning = readParValuesFromFitDir(bins, binningVars)  # append parameter values
        if parNames is not None:
          assert parNamesInBinning == parNames, f"The parameter set {parNamesInBinning} of dataset '{dataSet}' and binning '{binningVars}' is different from the parameter set {parNames} of the previous one"
        else:
          parNames = parNamesInBinning
    if binVarNames is not None:
      assert binVarNamesInDataSet == binVarNames, f"The binning variables {binVarNamesInDataSet} of dataset '{dataSet}' are different from the binning variables {binVarNames} of the previous one"
    else:
      binVarNames = binVarNamesInDataSet

  # plot fit parameters as 1D function of binning variable
  makePlots.setupPlotStyle()
  if parNames is not None:
    for parName in parNames:
      for binningVars in binVarNames:
        if len(binningVars) == 1:
          plotParValue1D(parValues, parName, binningVars, outputDirName)
