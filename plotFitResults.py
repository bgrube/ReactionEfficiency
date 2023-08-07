#!/usr/bin/env python3


import argparse
import collections
import ctypes
from dataclasses import dataclass  # builtin in Python 3.7+
import functools
import itertools
import numpy as np
import os
import sys
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

from uncertainties import UFloat, ufloat

import ROOT
if __name__ == "__main__":
  ROOT.PyConfig.DisableRootLogon = True  # do not change style of canvases loaded from fit result files  # type: ignore

import makePlots


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


BINNING_VAR_PLOT_INFO: Dict[str, Dict[str, str]] = {
  "BeamEnergy"         : {"label" : "E_{beam}",                      "unit" : "GeV"},
  "MissingProtonP"     : {"label" : "#it{p}^{miss}_{kin. fit}",      "unit" : "GeV/c"},
  "MissingProtonTheta" : {"label" : "#it{#theta}^{miss}_{kin. fit}", "unit" : "deg"},
  "MissingProtonPhi"   : {"label" : "#it{#phi}^{miss}_{kin. fit}",   "unit" : "deg"},
}

REMOVE_PARAM_BOX = False


@dataclass
class BinInfo:
  '''Stores information about a single kinematic bin'''
  name:    str               # bin name
  centers: Dict[str, float]  # dict with bin centers { <binning var> : <bin center>, ..., }
  dirName: str               # directory name for bins

  @property
  def varNames(self) -> Tuple[str, ...]:
    '''Returns names of kinematic variables used for binning'''
    return tuple(sorted(self.centers.keys()))

  @property
  def fitResultFileName(self) -> str:
    '''Returns fit result file name for this bin'''
    return f"{self.dirName}/ResultsHSMinuit2.root"

  def isSameBinAs(self, other: 'BinInfo') -> bool:
    '''Returns whether 2 BinInfo objects represent the same kinematic bin'''
    return (self.name == other.name) and (self.centers == other.centers)


@dataclass
class BinningInfo:
  '''Stores information about one particular kinematic binning'''
  infos:   List[BinInfo]  # info for all bins of the kinematic binning
  dirName: str            # directory that contains all bins of the kinematic binning

  @property
  def names(self) -> Tuple[str, ...]:
    '''Returns tuple with bin names'''
    return tuple(binInfo.name for binInfo in self.infos)

  @property
  def dirNames(self) -> Tuple[str, ...]:
    '''Returns tuple with directory names for all bins'''
    return tuple(binInfo.dirName for binInfo in self.infos)

  @property
  def varNames(self) -> Tuple[str, ...]:
    '''Returns names of kinematic variables used for binning'''
    varNames = None
    for binInfo in self.infos:
      varNamesInBin = binInfo.varNames
      if varNames is None:
        varNames = varNamesInBin
      else:
        assert varNamesInBin == varNames, f"Bins have inconsistent set of binning variables: bin '{binInfo.name}': {varNamesInBin} vs. {varNames} in previous bin"
    return varNames or tuple()

  @property
  def fitResultFileNames(self) -> Tuple[str, ...]:
    '''Returns tuple with fit-result file names for all bin'''
    return tuple(binInfo.fitResultFileName for binInfo in self.infos)


def getBinningFromDir(fitResultDirName: str) -> Union[BinningInfo, None]:
  '''Reads binning info from given directory'''
  binningFileName = f"{fitResultDirName}/DataBinsConfig.root"
  if not os.path.isfile(binningFileName):
    return None
  print(f"Loading binning from file '{binningFileName}'")
  bins = ROOT.Bins("HSBins", binningFileName)  # type: ignore
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


@dataclass
class ParInfo:
  '''Stores information about parameter values in a single kinematic bin'''
  binInfo: BinInfo            # info for the kinematic bin
  values:  Dict[str, UFloat]  # mapping of parameter names to values

  @property
  def names(self) -> Tuple[str, ...]:
    '''Returns parameter names'''
    return tuple(sorted(self.values.keys()))


def readParInfoForBin(
  binInfo:           BinInfo,
  fitParNamesToRead: Optional[Mapping[str, str]] = None,  # if dict { <new par name> : <par name>, ... } is set, only the given parameters are read, where `new par name` is the key used in the output
) -> ParInfo:
  '''Reads parameter values from fit result for given kinematic bin; an optional mapping selects parameters to read and translates parameter names'''
  fitResultFileName = binInfo.fitResultFileName
  print(f"Reading fit result object 'MinuitResult' from file '{fitResultFileName}'")
  fitResultFile  = ROOT.TFile.Open(fitResultFileName, "READ")  # type: ignore
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
  print("    Fit result status:", ", ".join([f"[{i}] {fitResult.statusCodeHistory(i)} = {fitResult.statusLabelHistory(i)}" for i in range(fitResult.numStatusHistory())]))
  fitResultFile.Close()
  return ParInfo(binInfo, parValuesInBin)


def readParInfosForBinning(binningInfo: BinningInfo) -> List[ParInfo]:
  '''Reads parameter values from fit results for given kinematic binning'''
  parInfos = []
  parNames = None  # used to compare parameter names for current and previous bin
  for binInfo in binningInfo.infos:
    parInfo = readParInfoForBin(binInfo)
    # ensure that the parameter sets of the fit functions are the same in all kinematic bins
    parNamesInBin = parInfo.names
    if parNames is not None:
      assert parNamesInBin == parNames, f"The parameter set {parNamesInBin} of this bin is different from the parameter set {parNames} of the previous one"
    parNames = parNamesInBin
    print(f"Read parameter values for kinematic bin: {parInfo}")
    parInfos.append(parInfo)
  return parInfos


def drawZeroLine(obj, style = ROOT.kDashed, color = ROOT.kBlack) -> None:  # type: ignore
  '''Helper function that draws zero line when necessary'''
  objType = obj.IsA().GetName()
  if (objType == "TCanvas") or (objType == "TPad"):
    xMin = ctypes.c_double()
    xMax = ctypes.c_double()
    yMin = ctypes.c_double()
    yMax = ctypes.c_double()
    obj.GetRangeAxis(xMin, yMin, xMax, yMax)
    if (yMin.value < 0) and (yMax.value > 0):
      zeroLine = ROOT.TLine()  # type: ignore
      zeroLine.SetLineStyle(style)
      zeroLine.SetLineColor(color)
      return zeroLine.DrawLine(xMin, 0, xMax, 0)
  elif objType.startswith("TH") or objType.startswith("TGraph") or objType.startswith("TMulti"):
    xAxis = obj.GetXaxis()
    yAxis = obj.GetYaxis()
    if (yAxis.GetXmin() < 0) and (yAxis.GetXmax() > 0):
      zeroLine = ROOT.TLine()  # type: ignore
      zeroLine.SetLineStyle(style)
      zeroLine.SetLineColor(color)
      return zeroLine.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 0, xAxis.GetBinUpEdge(xAxis.GetLast()), 0)
    elif (yAxis.GetXmin() > 0) and (yAxis.GetXmax() > 0):
      obj.SetMinimum(0)
  else:
    raise TypeError(f"drawZeroLine() not (yet) implemented for object of type '{objType}'")


def plotFitResult(
  binInfo:     BinInfo,
  fitVariable: str,
  pdfDirName:  Optional[str] = None,  # overrides default PDF output path (i.e. same dir as fit result file) if set
) -> None:
  '''Plots fit result for given kinematic bin'''
  fitResultFileName = binInfo.fitResultFileName
  if not os.path.isfile(fitResultFileName):
    print(f"Cannot find file '{fitResultFileName}'; skipping")
    return
  print(f"Plotting fit result in file '{fitResultFileName}'")
  fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")  # type: ignore
  canvName = f"{binInfo.name or ''}_{fitVariable}"
  canv = fitResultFile.Get(canvName)
  # improve TPaveText with fit parameters
  dataFitPad = canv.GetListOfPrimitives().FindObject(f"{canvName}_1")
  paramBox = dataFitPad.GetListOfPrimitives().FindObject(f"{binInfo.name or ''}TotalPDF_paramBox")
  if REMOVE_PARAM_BOX:
    # remove box completely
    dataFitPad.GetListOfPrimitives().Remove(paramBox)
  else:
    # only remove filled frame
    paramBox.SetBorderSize(0)
    paramBox.SetFillStyle(0)
  drawZeroLine(dataFitPad)
  pdfFileName = ("Overall" if binInfo.name == "" else "") + f"{canv.GetName()}.pdf"
  if pdfDirName:
    canv.SaveAs(f"{pdfDirName}/{pdfFileName}")
  else:
    canv.SaveAs(f"{binInfo.dirName}/{pdfFileName}")
  fitResultFile.Close()


def plotFitResults(
  binningInfo: BinningInfo,
  fitVariable: str,
) -> None:
  '''Plots fit results for all kinematic bins'''
  for binInfo in binningInfo.infos:
    plotFitResult(binInfo, fitVariable, pdfDirName = binningInfo.dirName)


def getParValuesForGraph1D(
  binVarName: str,  # name of the binning variable, i.e. x-axis
  parName:    str,  # name of value, i.e. y-axis
  parInfos:   Sequence[ParInfo],
) -> Tuple[Tuple[float, UFloat]]:
  '''Extracts information needed to plot parameter with given name as a function of the given bin variable from list of ParInfos'''
  graphValues: Tuple[Tuple[float, UFloat]] = tuple((parInfo.binInfo.centers[binVarName], parInfo.values[parName])
    for parInfo in parInfos if (binVarName in parInfo.binInfo.varNames) and (parName in parInfo.names))
  return graphValues


def getParValueGraph1D(
  graphValues:     Sequence[Tuple[float, UFloat]],
  shiftByFraction: float = 0,
) -> Any:  #TODO there does not seem to be a way to specify ROOT types
  '''Creates ROOT.TGraphErrors from given values'''
  if not graphValues:
    print("No data to plot")
    return
  xVals = np.array([graphVal[0]               for graphVal in graphValues], dtype = "d")
  yVals = np.array([graphVal[1].nominal_value for graphVal in graphValues], dtype = "d")
  yErrs = np.array([graphVal[1].std_dev       for graphVal in graphValues], dtype = "d")
  # shift x values by fraction of total x range
  if shiftByFraction != 0:
    xRange = max(xVals) - min(xVals)
    shift  = xRange * shiftByFraction
    xVals = xVals + shift
  # report weighted average
  meanEff = np.average(yVals, weights = [1 / (yErr**2) for yErr in yErrs])
  print(f"    weighted mean of efficiencies = {meanEff}")
  return ROOT.TGraphErrors(len(xVals), xVals, yVals, ROOT.nullptr, yErrs)  # type: ignore


def plotParValue1D(
  parInfos:          Mapping[str, Sequence[ParInfo]],
  parName:           str,  # name of parameter to plot
  binningVar:        str,  # name of binning variable to plot
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNameSuffix: str   = "",
  particle:          str   = "Proton",
  channel:           str   = "4pi",
  markerSize:        float = 0.75,
) -> None:
  '''Overlay parameter values for datasets for 1-dimensional binning for given parameter for given binning variable'''
  print(f"Plotting parameter '{parName}' as a function of binning variable '{binningVar}'")
  parValueMultiGraph = ROOT.TMultiGraph()  # type: ignore
  parValueGraphs = {}  # store graphs here to keep them in memory
  shiftFraction = 0
  styleIndex = 0
  for dataSet in parInfos:
    graph = parValueGraphs[dataSet] = getParValueGraph1D(getParValuesForGraph1D(binningVar, parName, parInfos[dataSet]), shiftFraction)
    shiftFraction += 0.01
    graph.SetTitle(dataSet)  # type: ignore
    makePlots.setCbFriendlyStyle(graph, styleIndex, skipBlack = False)
    styleIndex += 1
    parValueMultiGraph.Add(graph)
  parValueMultiGraph.SetTitle(f"Fit parameter {parName}, {particle} ({channel})")
  assert binningVar in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVar}'"
  parValueMultiGraph.GetXaxis().SetTitle(f"{BINNING_VAR_PLOT_INFO[binningVar]['label']} ({BINNING_VAR_PLOT_INFO[binningVar]['unit']})")
  parValueMultiGraph.GetYaxis().SetTitle(parName)
  canv = ROOT.TCanvas(f"{particle}_{channel}_{parName}_{binningVar}{pdfFileNameSuffix}", "")  # type: ignore
  parValueMultiGraph.Draw("APZ")
  hist = parValueMultiGraph.GetHistogram()
  plotRange = hist.GetMaximum() - hist.GetMinimum()
  minVal = min(tuple(ROOT.TMath.MinElement(graph.GetN(), graph.GetY())  for graph in parValueMultiGraph.GetListOfGraphs()))  # type: ignore
  maxVal = max(tuple(ROOT.TMath.MaxElement(graph.GetN(), graph.GetY())  for graph in parValueMultiGraph.GetListOfGraphs()))  # type: ignore
  maxErr = max(tuple(ROOT.TMath.MaxElement(graph.GetN(), graph.GetEY()) for graph in parValueMultiGraph.GetListOfGraphs()))  # type: ignore
  if 2 * maxErr / plotRange > 0.9:  # one error bar dominates plot range
    plotRange = 2 * (maxVal - minVal)  # new plot range = 2 * value range
    plotRangeCenter = (maxVal + minVal) / 2
    minVal = plotRangeCenter - plotRange / 2
    maxVal = plotRangeCenter + plotRange / 2
    print(f"Adjusting plot range to [{minVal}, {maxVal}]")
    parValueMultiGraph.SetMinimum(minVal)
    parValueMultiGraph.SetMaximum(maxVal)
    canv.Modified()
    canv.Update()
  legend = canv.BuildLegend()
  legend.SetFillStyle(0)
  legend.SetBorderSize(0)
  drawZeroLine(parValueMultiGraph)
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


if __name__ == "__main__":
  makePlots.printGitInfo()
  ROOT.gROOT.SetBatch(True)  # type: ignore
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")  # type: ignore

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Plots BruFit results in given directory.")
  parser.add_argument("outputDirName", type = str, nargs = "?", default = "./BruFitOutput", help = "The path to the BruFit output directory; (default: '%(default)s')")
  args = parser.parse_args()
  dataSets    = ["Total", "Found", "Missing"]
  fitVariable = "MissingMassSquared_Measured"

  parInfos:    Dict[str, List[ParInfo]]           = {}    # parInfos[<dataset>][<bin>]
  parNames:    Union[Tuple[str, ...], None]       = None
  binVarNames: Union[List[Tuple[str, ...]], None] = None  # binning variables for each binning
  for dataSet in dataSets:
    fitResultDirName  = f"{args.outputDirName}/{dataSet}"
    print(f"Plotting overall fit result for '{dataSet}' dataset")
    plotFitResult(BinInfo("", {}, fitResultDirName), fitVariable)
    binVarNamesInDataSet = []
    for binningInfo in getBinningInfosFromDir(fitResultDirName):
      if binningInfo:
        binVarNamesInDataSet.append(binningInfo.varNames)
        if binningInfo.dirNames:
          print(f"Plotting fit results for binning variable(s) '{binningInfo.varNames}' for '{dataSet}' dataset")
          plotFitResults(binningInfo, fitVariable)
          if not dataSet in parInfos:
            parInfos[dataSet] = []
          parInfos[dataSet][len(parInfos[dataSet]):] = readParInfosForBinning(binningInfo)  # append parameter values
          parNamesInBinning = parInfos[dataSet][-1].names  # readParInfosForBinning() ensures that parameter names are identical within a binning
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
          plotParValue1D(parInfos, parName, binningVars[0], args.outputDirName)
