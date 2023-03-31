#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT

import argparse
import array
import collections
import ctypes
from dataclasses import dataclass  # builtin in Python 3.7+
import itertools
import os
import sys
from typing import Optional, Union, Dict, Tuple, List, Iterable, Sequence, Mapping

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


@dataclass
class BinInfo:
  '''Stores information about a single kinematic bin'''
  name:    Optional[str]              = None  # bin name
  centers: Optional[Dict[str, float]] = None  # dict with bin centers { <binning var> : <bin center>, ..., }
  dirName: Optional[str]              = None  # directory names for all bins ( <bin dir name>, ... )

  @property
  def varNames(self) -> Union[Tuple[str, ...], None]:
    '''Returns names of kinematic variables used for binning'''
    return None if self.centers is None else tuple(sorted(self.centers.keys()))

  @property
  def fitResultFileName(self) -> Union[str, None]:
    '''Returns fit result file name for this bin'''
    return None if self.dirName is None else f"{self.dirName}/ResultsHSMinuit2.root"


@dataclass
class BinningInfo:
  '''Stores information about one particular kinematic binning'''
  binInfos: Optional[Tuple[BinInfo, ...]] = None
  dirName:  Optional[str]                 = None

  @property
  def names(self) -> Union[Tuple[Union[str, None], ...], None]:
    '''Returns tuple with bin names'''
    return None if self.binInfos is None else tuple(binInfo.name for binInfo in self.binInfos)

  @property
  def dirNames(self) -> Union[Tuple[Union[str, None], ...], None]:
    '''Returns tuple with directory names for all bins'''
    return None if self.binInfos is None else tuple(binInfo.dirName for binInfo in self.binInfos)

  @property
  def varNames(self) -> Union[Tuple[str, ...], None]:
    '''Returns names of kinematic variables used for binning'''
    varNames = None
    if self.binInfos:
      for binInfo in self.binInfos:
        varNamesInBin = binInfo.varNames
        if varNames is None:
          varNames = varNamesInBin
        else:
          assert varNamesInBin == varNames, f"Bins have inconsistent set of binning variables: bin '{binInfo.name}': {varNamesInBin} vs. {varNames} in previous bin"
    return varNames

  @property
  def fitResultFileNames(self) -> Union[Tuple[Union[str, None]], None]:
    '''Returns tuple with fit-result file names for all bin'''
    return None if self.binInfos is None else tuple(binInfo.fitResultFileName for binInfo in self.binInfos)


def getBinningFromDir(fitResultDirName: str) -> Union[BinningInfo, None]:
  '''Reads binning info from given directory'''
  binningFileName = f"{fitResultDirName}/DataBinsConfig.root"
  if not os.path.isfile(binningFileName):
    return None
  print(f"Loading binning from file '{binningFileName}'")
  bins = ROOT.Bins("HSBins", binningFileName)
  print("Found binning:")
  bins.PrintAxis()
  binNames = tuple(str(binName) for binName in bins.GetBinNames())
  axes = bins.GetVarAxis()
  binVarNames = tuple(axis.GetName() for axis in axes)
  axisBinIndexRanges = tuple(range(1, axis.GetNbins() + 1) for axis in axes)
  binInfos: List[BinInfo] = []
  for axisBinIndices in itertools.product(*axisBinIndexRanges):  # loop over all tuples of bin indices for the axes
    axisBinCenters = tuple(axes[axisIndex].GetBinCenter(axisBinIndex) for axisIndex, axisBinIndex in enumerate(axisBinIndices))
    binIndex = bins.FindBin(*axisBinCenters)  #!Note! the unpacking works only for up to 6 binning dimensions
    binName = binNames[binIndex]
    binInfos.append(BinInfo(
      name = binName,
      centers = {binVarNames[axisIndex] : axisBinCenter for axisIndex, axisBinCenter in enumerate(axisBinCenters)},
      dirName = f"{fitResultDirName}/{str(binName)}",
    ))
  return BinningInfo(binInfos, fitResultDirName)


def getBinningInfosFromDir(fitResultDirName: str) -> List[Union[BinningInfo, None]]:
  '''Reads binning infos from all 1st-level subdirectories that contain a 'DataBinsConfig.root' file'''
  binningInfos: List[Union[BinningInfo, None]] = []
  # find all subdirectories with binning files non-recursively
  subDirNames: List[str] = sorted([entry.path for entry in os.scandir(fitResultDirName) if entry.is_dir() and os.path.isfile(f"{entry.path}/DataBinsConfig.root")])
  # get binning info from all subdirectories
  for subDirName in subDirNames:
    print(f"Found binning info in directory '{subDirName}'")
    binningInfos.append(getBinningFromDir(subDirName))
  return binningInfos


#TODO convert to dataclass
#TODO complete docstrings
#TODO complete type annotations
# reads fit result from given file and returns
#    dict with parameter values { <par name> : <par value>, ... }
#    tuple with parameter names (<par name>, ... )
def readParValuesFromFitFile(
  fitResultFileName: Union[str, None],
  fitParNamesToRead: Optional[Mapping[str, str]] = None  # if dict { <new par name> : <par name>, ... } is set, only the given parameters are read, where `new par name` is the key used in the output
):
  if not fitResultFileName:
    return
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


# reads parameter values from fit results for given binning
#     list of dicts with parameter values [ { <binning var> : <bin center>, ..., <par name> : <par value>, ... }, ... ]
#     tuple with parameter names (<par name>, ... )
def readParValuesForBinning(binningInfo: BinningInfo):
  parValues = []
  parNames = None  # used to compare parameter names in current and previous bin
  if binningInfo.binInfos:
    for binInfo in binningInfo.binInfos:
      parValuesInBin, parNamesInBin = readParValuesFromFitFile(binInfo.fitResultFileName)
      # ensure that the parameter sets of the fit functions are the same in all kinematic bins
      if parNames is not None:
        assert parNamesInBin == parNames, f"The parameter set {parNamesInBin} of this bin is different from the parameter set {parNames} of the previous one"
      parNames = parNamesInBin
      # copy axis name(s) and bin center(s)
      parValuesInBin.update(binInfo.centers)
      print(f"Read parameter values for kinematic bin: {parValuesInBin}")
      parValues.append(parValuesInBin)
  return parValues, parNames


def drawZeroLine(obj, style = ROOT.kDashed, color = ROOT.kBlack) -> None:
  '''Helper function that draws zero line when necessary'''
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


def plotFitResult(
  fitResultDirName: Union[str, None],
  fitVariable:      str,
  binName:          Union[str, None] = "",
  pdfDirName:       Optional[str]    = None  # overrides default PDF output path (i.e. same dir as fit result file) if set
) -> None:
  '''Plots fit result in given directory'''
  if not fitResultDirName or not binName:
    return
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
  if pdfDirName:
    canv.SaveAs(f"{pdfDirName}/{pdfFileName}")
  else:
    canv.SaveAs(f"{fitResultDirName}/{pdfFileName}")
  fitResultFile.Close()


def plotFitResults(
  binNames:          Union[Sequence[Union[str, None]], None],
  fitResultDirNames: Union[Sequence[Union[str, None]], None],
  fitVariable:       str,
  pdfDirName:        Optional[str] = None
):
  '''Plots fit results for all kinematic bins'''
  if fitResultDirNames and binNames:
    for index, fitResultDirName  in enumerate(fitResultDirNames):
      plotFitResult(fitResultDirName, fitVariable, binNames[index], pdfDirName)


# returns
#     arrays of x, y, and y-uncertainty values
#     x-axis label and unit
def getDataPointArrays1D(
  binVarName: str,  # name of the binning variable, i.e. x-axis
  valueName:  str,  # name of value, i.e. y-axis
  values            # list of dicts with data values [ { <binning var> : <bin center>, ..., <value name> : <value> }, ... ]
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


def plotParValue1D(
  parValues,      # dict of lists of dicts with parameter values { <dataset> : [ { <binning var> : <bin center>, ..., <par name> : <par value>, ... }, ... ], ... }
  parName,        # name of parameter to plot
  binningVars,    # tuple with binning variables (<binning var>, ... )
  pdfDirName,     # directory name the PDF file will be written to
  pdfFileNameSuffix = "",
  particle          = "Proton",
  channel           = "4pi",
  markerSize        = 0.75
):
  '''Overlay parameter values for datasets for 1-dimensional binning'''
  binningVar  = binningVars[0]
  print(f"Plotting parameter '{parName}' as a function of binning variable '{binningVar}'")
  parValueMultiGraph = ROOT.TMultiGraph()
  parValueGraphs = {}  # store graphs here to keep them in memory
  for dataSet in parValues:
    xVals, yVals, yErrs, binVarLabel, binVarUnit = getDataPointArrays1D(binningVar, parName, parValues[dataSet])
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
  canv = ROOT.TCanvas(f"{particle}_{channel}_{parName}_{binningVar}{pdfFileNameSuffix}", "")
  parValueMultiGraph.Draw("APZ")
  canv.BuildLegend()
  drawZeroLine(parValueMultiGraph)
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Plots BruFit results.")
  parser.add_argument("outputDirName", type = str, nargs = "?", default = "./BruFitOutput", help = "The path to the BruFit output directory; (default: '%(default)s')")
  args = parser.parse_args()
  dataSets    = ["Total", "Found", "Missing"]
  fitVariable = "MissingMassSquared_Measured"  #TODO get this info from BruFit's ROOT.Setup

  parValues   = {}  # dict of lists of dicts with parameter values { <dataset> : [ { <binning var> : <bin center>, ..., <par name> : <par value>, ... }, ... ], ... }
  parNames    = None  # tuple with parameter names (<par name>, ... )
  binVarNames = None  # list of lists with binning variables for each binning [ [ <binning var>, ... ], ... ]
  for dataSet in dataSets:
    fitResultDirName  = f"{args.outputDirName}/{dataSet}"
    print(f"Plotting overall fit result for '{dataSet}' dataset")
    plotFitResult(fitResultDirName, fitVariable)
    binVarNamesInDataSet = []
    for binningInfo in getBinningInfosFromDir(fitResultDirName):
      if binningInfo:
        binVarNamesInDataSet.append(binningInfo.varNames)
        if binningInfo.dirNames:
          print(f"Plotting fit results for binning variable(s) '{binningInfo.varNames}' for '{dataSet}' dataset")
          plotFitResults(
            binNames          = binningInfo.names,
            fitResultDirNames = binningInfo.dirNames,
            fitVariable       = fitVariable,
            pdfDirName        = binningInfo.dirName,
          )
          if not dataSet in parValues:
            parValues[dataSet] = []
          parValues[dataSet][len(parValues[dataSet]):], parNamesInBinning = readParValuesForBinning(binningInfo)  # append parameter values
          if parNames is not None:
            assert parNamesInBinning == parNames, f"The parameter set {parNamesInBinning} of '{dataSet}' dataset and binning '{binningInfo.varNames}' is different from the parameter set {parNames} of the previous one"
          else:
            parNames = parNamesInBinning
    if binVarNames is not None:
      assert binVarNamesInDataSet == binVarNames, f"The binning variables {binVarNamesInDataSet} of '{dataSet}' dataset are different from the binning variables {binVarNames} of the previous one"
    else:
      binVarNames = binVarNamesInDataSet

  # plot fit parameters as 1D function of binning variable
  makePlots.setupPlotStyle()
  if parNames and binVarNames:
    for parName in parNames:
      for binningVars in binVarNames:
        if len(binningVars) == 1:
          plotParValue1D(parValues, parName, binningVars, args.outputDirName)
