#!/usr/bin/env python3


import argparse
import collections.abc
from dataclasses import dataclass
import functools
import itertools
import math
import os
import sys
from typing import (
  Any,
  Dict,
  List,
  Mapping,
  Optional,
  Sequence,
  Tuple,
  Union,
)  #TODO use A | B syntax with `from __future__ import annotations` for Union and Optional, see https://adamj.eu/tech/2022/10/17/python-type-hints-old-and-new-syntaxes/; use builtin types, see https://peps.python.org/pep-0585/#implementation

from uncertainties import UFloat, ufloat

import ROOT
if __name__ == "__main__":
  ROOT.PyConfig.DisableRootLogon = True  # do not change style of canvases loaded from fit result files

# import plotBeautifiers
from plotTools import (
  drawZeroLine,
  getGraph1DFromValues,
  getGraph2DFromValues,
  printGitInfo,
  setCbFriendlyStyle,
  setupPlotStyle,
)


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


BINNING_VAR_PLOT_INFO: Dict[str, Dict[str, str]] = {
  "BeamEnergy"         : {"label" : "#it{E}_{beam}",                 "unit" : "GeV"},
  "MissingProtonP"     : {"label" : "#it{p}_{miss}^{kin. fit}",      "unit" : "GeV/#it{c}"},
  "MissingProtonTheta" : {"label" : "#it{#theta}_{miss}^{kin. fit}", "unit" : "deg"},
  "MissingProtonPhi"   : {"label" : "#it{#phi}_{miss}^{kin. fit}",   "unit" : "deg"},
}

def getAxisInfoForBinningVar(binningVar: Union[str, Tuple[str, str]]) -> Tuple[str, str, str]:
  varName = None
  if isinstance(binningVar, str):
    varName = binningVar
  elif isinstance(binningVar, Sequence):
    varName = binningVar[0]
  assert varName is not None and varName in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{varName}'"
  varLabel = BINNING_VAR_PLOT_INFO[varName]['label']
  varUnit  = BINNING_VAR_PLOT_INFO[varName]['unit']
  return (varName, varLabel, varUnit)


REMOVE_PARAM_BOX = False


@dataclass
class BinInfo:
  """Stores information about a single kinematic bin"""
  name:    str               # bin name
  centers: Dict[str, float]  # bin centers { <binning var> : <bin center>, ... }
  widths:  Dict[str, float]  # bin widths  { <binning var> : <bin width>,  ... }
  dirName: str               # directory name for bins

  @property
  def varNames(self) -> Tuple[str, ...]:
    """Returns names of kinematic variables used for binning"""
    return tuple(sorted(self.centers.keys()))

  @property
  def fitResultFileName(self) -> str:
    """Returns fit result file name for this bin"""
    return f"{self.dirName}/ResultsHSMinuit2.root"

  def isSameBinAs(self, other: 'BinInfo') -> bool:
    """Returns whether 2 BinInfo objects represent the same kinematic bin"""
    return (self.name == other.name) and (self.centers == other.centers)


@dataclass
class BinningInfo:
  """Stores information about one particular kinematic binning"""
  infos:   List[BinInfo]  # info for all bins of the kinematic binning
  dirName: str            # directory that contains all bins of the kinematic binning

  @property
  def names(self) -> Tuple[str, ...]:
    """Returns tuple with bin names"""
    return tuple(binInfo.name for binInfo in self.infos)

  @property
  def dirNames(self) -> Tuple[str, ...]:
    """Returns tuple with directory names for all bins"""
    return tuple(binInfo.dirName for binInfo in self.infos)

  @property
  def varNames(self) -> Tuple[str, ...]:
    """Returns names of kinematic variables used for binning"""
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
    """Returns tuple with fit-result file names for all bin"""
    return tuple(binInfo.fitResultFileName for binInfo in self.infos)


def getBinningFromDir(fitResultDirName: str) -> Optional[BinningInfo]:
  """Reads binning info from given directory"""
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
    axisBinWidths  = tuple(axes[axisIndex].GetBinWidth (axisBinIndex) for axisIndex, axisBinIndex in enumerate(axisBinIndices))
    binIndex = bins.FindBin(*axisBinCenters)  #!Note! the unpacking works only for up to 6 binning dimensions
    binName = binNames[binIndex]
    binInfos.append(BinInfo(
      name    = binName,
      centers = {binVarNames[axisIndex] : axisBinCenter for axisIndex, axisBinCenter in enumerate(axisBinCenters)},
      widths  = {binVarNames[axisIndex] : axisBinWidth  for axisIndex, axisBinWidth  in enumerate(axisBinWidths )},
      dirName = f"{fitResultDirName}/{str(binName)}",
    ))
  return BinningInfo(binInfos, fitResultDirName)


def getBinningInfosFromDir(fitResultDirName: str) -> List[Optional[BinningInfo]]:
  """Reads binning infos from all 1st-level subdirectories that contain a 'DataBinsConfig.root' file"""
  binningInfos: List[Optional[BinningInfo]] = []
  # find all subdirectories with binning files non-recursively
  subDirNames: List[str] = sorted([entry.path for entry in os.scandir(fitResultDirName) if entry.is_dir() and os.path.isfile(f"{entry.path}/DataBinsConfig.root")])
  # get binning info from all subdirectories
  for subDirName in subDirNames:
    print(f"Found binning info in directory '{subDirName}'")
    binningInfos.append(getBinningFromDir(subDirName))
  return binningInfos


@dataclass
class ParInfo:
  """Stores information about parameter values in a single kinematic bin"""
  binInfo: BinInfo            # info for the kinematic bin
  values:  Dict[str, UFloat]  # mapping of parameter names to values

  @property
  def names(self) -> Tuple[str, ...]:
    """Returns parameter names"""
    return tuple(sorted(self.values.keys()))


def readParInfoForBin(
  binInfo:           BinInfo,
  fitParNamesToRead: Optional[Mapping[str, str]] = None,  # if dict { <new par name> : <par name>, ... } is set, only the given parameters are read, where `new par name` is the key used in the output
) -> ParInfo:
  """Reads parameter values from fit result for given kinematic bin; an optional mapping selects parameters to read and translates parameter names"""
  fitResultFileName = binInfo.fitResultFileName
  if not os.path.isfile(fitResultFileName):
    print(f"Cannot find file '{fitResultFileName}'. Skipping bin {binInfo}.")
    return ParInfo(binInfo, {})
  print(f"Reading fit result object 'MinuitResult' from file '{fitResultFileName}'")
  fitResultFile  = ROOT.TFile.Open(fitResultFileName, "READ")
  fitResult      = fitResultFile.Get("MinuitResult")
  fitPars        = fitResult.floatParsFinal()
  parValuesInBin = {}
  if isinstance(fitParNamesToRead, collections.abc.Mapping):
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
  """Reads parameter values from fit results for given kinematic binning"""
  parInfos: List[ParInfo] = []
  parNames = None  # used to compare parameter names for current and previous bin
  for binInfo in binningInfo.infos:
    parInfo = readParInfoForBin(binInfo)
    if parInfo.values:
      # ensure that the parameter sets of the fit functions are the same in all kinematic bins, for which parameters could be read successfully
      parNamesInBin = parInfo.names
      if parNames is not None:
        assert parNamesInBin == parNames, f"The parameter set {parNamesInBin} of this bin is different from the parameter set {parNames} of the previous one"
      parNames = parNamesInBin
      print(f"Read parameter values for kinematic bin: {parInfo}")
      parInfos.append(parInfo)
  return parInfos


def plotFitResult(
  binInfo:     BinInfo,
  fitVariable: str,
  pdfDirName:  Optional[str] = None,  # overrides default PDF output path (i.e. same dir as fit result file) if set
) -> None:
  """Plots fit result for given kinematic bin"""
  fitResultFileName = binInfo.fitResultFileName
  if not os.path.isfile(fitResultFileName):
    print(f"Cannot find file '{fitResultFileName}'. Skipping bin {binInfo}.")
    return
  print(f"Plotting fit result in file '{fitResultFileName}'")
  fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
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
  """Plots fit results for all kinematic bins"""
  for binInfo in binningInfo.infos:
    plotFitResult(binInfo, fitVariable, pdfDirName = binningInfo.dirName)


def getParValuesForGraph1D(
  binVarName: str,  # name of the binning variable, i.e. x-axis
  parName:    str,  # name of value, i.e. y-axis
  parInfos:   Sequence[ParInfo],
) -> Tuple[Tuple[UFloat, UFloat], ...]:
  """Extracts information needed to plot parameter with given name as a function of the given bin variable from list of ParInfos"""
  graphValues: Tuple[Tuple[UFloat, UFloat], ...] = tuple(
    (ufloat(parInfo.binInfo.centers[binVarName], parInfo.binInfo.widths[binVarName] / 2.0), parInfo.values[parName])
    for parInfo in parInfos
    if (binVarName in parInfo.binInfo.varNames) and (parName in parInfo.names) and (len(parInfo.binInfo.varNames) == 1)
  )
  return graphValues


#TODO add sequence of beautifier functors as optional argument
def plotGraphs1D(
  graphOrGraphs:     Union[ROOT.TGraph, Sequence[Tuple[str, ROOT.TGraph]]],
  binningVar:        str,
  yAxisTitle:        str,
  pdfDirName:        str,
  pdfFileBaseName:   str,
  pdfFileNamePrefix: str = "Proton_4pi_",
  pdfFileNameSuffix: str = "",
  graphTitle:        Optional[str] = None,
  graphMinimum:      Optional[float] = 0.0,
  graphMaximum:      Optional[float] = None,
  skipBlack:         bool = True,
  drawLegend:        Optional[bool] = None,
  forceXRange:       Optional[Tuple[float, float]] = None,
  beautifiers:       Sequence[Any] = [],  #TODO improve type hint by making a base class for beautifiers
) -> None:
  """Generic function that plots the given graph(s)"""
  graphs: List[Tuple[str, ROOT.TGraph]] = []
  if isinstance(graphOrGraphs, ROOT.TGraph):
    graphs.append(("", graphOrGraphs))
  elif isinstance(graphOrGraphs, Sequence):
    graphs = list(graphOrGraphs)
  else:
    raise TypeError(f"`graphOrGraphs` must be an instance of either ROOT.TGraph or Sequence but is a {type(graphOrGraphs)}")
  multiGraph = ROOT.TMultiGraph()
  for styleIndex, (legendLabel, graph) in enumerate(graphs):
    graph.SetTitle(legendLabel)
    setCbFriendlyStyle(graph, styleIndex, skipBlack = False if len(graphs) == 1 else skipBlack)
    multiGraph.Add(graph)
  if graphTitle is not None:
    multiGraph.SetTitle(graphTitle)  # !Note! if this is executed after setting axis titles, no title is printed; seems like a ROOT bug
  #TODO use getAxisInfoForBinningVar everywhere
  #TODO make this function more generic by moving this into a beautifier
  assert binningVar in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVar}'"
  multiGraph.GetXaxis().SetTitle(f"{BINNING_VAR_PLOT_INFO[binningVar]['label']} ({BINNING_VAR_PLOT_INFO[binningVar]['unit']})")
  multiGraph.GetYaxis().SetTitle(yAxisTitle)
  if graphMinimum is not None:
    multiGraph.SetMinimum(graphMinimum)
  if graphMaximum is not None:
    multiGraph.SetMaximum(graphMaximum)
  legendLabels = tuple(legendLabel.replace(' ', '_') for legendLabel, _ in graphs)  # reformat legend labels so that they can be used in PDF file name
  if len(legendLabels) == 1 and legendLabels[0] == "":  # only 1 graph and no legend label
    legendLabels = None
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{pdfFileBaseName}_{binningVar}"
                      + ("" if legendLabels is None else f"_{'_'.join(legendLabels)}") + pdfFileNameSuffix, "")
  multiGraph.Draw("APZ")
  if forceXRange is not None:
    canv.Modified()
    multiGraph.GetXaxis().SetLimits(forceXRange[0], forceXRange[1])
  for beautifier in beautifiers:
    beautifier.draw(multiGraph)
  if drawLegend == True or (drawLegend is None and legendLabels is not None):
    legend = ROOT.TLegend(0.3, 0.21, 0.3, 0.21)
    for graph in multiGraph.GetListOfGraphs():
      label = graph.GetTitle()
      if not label:
        label = graph.GetName()
      legend.AddEntry(graph, label, "LPF")
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.Draw()
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


def plotParValue1D(
  parInfos:          Mapping[str, Sequence[ParInfo]],
  parName:           str,  # name of parameter to plot
  binningVar:        str,  # name of binning variable to plot
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str = "Proton_4pi_",
  pdfFileNameSuffix: str = "",
) -> None:
  """Overlays values of given parameter for given datasets for 1-dimensional binning with given binning variables"""
  print(f"Plotting parameter '{parName}' as a function of binning variable '{binningVar}'")
  parValueGraphs: List[Tuple[str, ROOT.TGraph]] = [(dataSet, getGraph1DFromValues(getParValuesForGraph1D(binningVar, parName, parInfos[dataSet])))
                                                   for dataSet in parInfos]
  plotGraphs1D(
    graphOrGraphs     = parValueGraphs,
    binningVar        = binningVar,
    yAxisTitle        = parName,
    pdfDirName        = pdfDirName,
    pdfFileBaseName   = parName,
    pdfFileNamePrefix = pdfFileNamePrefix,
    pdfFileNameSuffix = pdfFileNameSuffix,
    graphTitle        = f"Fit parameter {parName})",
    skipBlack         = False,
  )
  #TODO needed?
  # hist = parValueMultiGraph.GetHistogram()
  # plotRange = hist.GetMaximum() - hist.GetMinimum()
  # minVal: float = min(tuple(ROOT.TMath.MinElement(graph.GetN(), graph.GetY())  for graph in parValueMultiGraph.GetListOfGraphs()))
  # maxVal: float = max(tuple(ROOT.TMath.MaxElement(graph.GetN(), graph.GetY())  for graph in parValueMultiGraph.GetListOfGraphs()))
  # maxErr: float = max(tuple(ROOT.TMath.MaxElement(graph.GetN(), graph.GetEY()) for graph in parValueMultiGraph.GetListOfGraphs()))
  # if 2 * maxErr / plotRange > 0.9:  # one error bar dominates plot range
  #   plotRange = 2 * (maxVal - minVal)  # new plot range = 2 * value range
  #   plotRangeCenter = (maxVal + minVal) / 2
  #   minVal = plotRangeCenter - plotRange / 2
  #   maxVal = plotRangeCenter + plotRange / 2
  #   print(f"Adjusting plot range to [{minVal}, {maxVal}]")
  #   parValueMultiGraph.SetMinimum(minVal)
  #   parValueMultiGraph.SetMaximum(maxVal)
  #   canv.Modified()
  #   canv.Update()
  #TODO pass to plotGraphs1D()
  # drawZeroLine(parValueMultiGraph)


def getParValuesForGraph2D(
  binVarNames: Sequence[str],  # names of the binning variables, i.e. x-axis and y-axis
  parName:     str,  # name of value, i.e. z-axis
  parInfos:    Sequence[ParInfo],
) -> Tuple[Tuple[UFloat, UFloat, UFloat], ...]:
  """Extracts information needed to plot parameter with given name as a function of the given bin variables from list of ParInfos"""
  graphValues: Tuple[Tuple[UFloat, UFloat, UFloat], ...] = tuple(
    (
      ufloat(parInfo.binInfo.centers[binVarNames[0]], parInfo.binInfo.widths[binVarNames[0]] / 2.0),
      ufloat(parInfo.binInfo.centers[binVarNames[1]], parInfo.binInfo.widths[binVarNames[1]] / 2.0),
      parInfo.values[parName]
    )
    for parInfo in parInfos
    if (binVarNames[0] in parInfo.binInfo.varNames) and (binVarNames[1] in parInfo.binInfo.varNames) and (parName in parInfo.names) and (len(parInfo.binInfo.varNames) == 2)
  )
  return graphValues


def plotParValue2D(
  parInfos:          Mapping[str, Sequence[ParInfo]],
  parName:           str,  # name of parameter to plot
  binningVars:       Sequence[str],  # names of binning variables to plot
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str = "Proton_4pi_",
  pdfFileNameSuffix: str = "",
) -> None:
  """Overlays values of given parameter for given datasets for 2-dimensional binning with given binning variables"""
  savedFrameFillColor = ROOT.gStyle.GetFrameFillColor()
  ROOT.gStyle.SetFrameFillColor(0)  # switch back to default; otherwise graphs obstruct histogram frame
  print(f"Plotting parameter '{parName}' as a function of binning variables '{binningVars}'")
  parValueGraphs: Dict[str, ROOT.TGraph2DErrors] = {}  # store graphs here to keep them in memory
  xMin = yMin = zMin = +math.inf
  xMax = yMax = zMax = -math.inf
  for dataSet in parInfos:
    graph = getGraph2DFromValues(getParValuesForGraph2D(binningVars, parName, parInfos[dataSet]))
    if graph is not None:
      parValueGraphs[dataSet] = graph
      xMin = min(xMin, graph.GetXminE())
      yMin = min(yMin, graph.GetYminE())
      zMin = min(zMin, graph.GetZmin())
      xMax = max(xMax, graph.GetXmaxE())
      yMax = max(yMax, graph.GetYmaxE())
      zMax = max(zMax, graph.GetZmax())
  if not parValueGraphs:  # nothing to plot
    return
  canv = ROOT.TCanvas()
  # construct dummy histogram that ensures same axes for all graphs
  assert binningVars[0] in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVars[0]}'"
  assert binningVars[1] in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVars[1]}'"
  hist = ROOT.TH3F(
    f"{pdfFileNamePrefix}{parName}_{binningVars[0]}_{binningVars[1]}{pdfFileNameSuffix}",
    (f";{BINNING_VAR_PLOT_INFO[binningVars[0]]['label']} ({BINNING_VAR_PLOT_INFO[binningVars[0]]['unit']})"
     f";{BINNING_VAR_PLOT_INFO[binningVars[1]]['label']} ({BINNING_VAR_PLOT_INFO[binningVars[1]]['unit']})"
     f";{parName}"),
    1, xMin, xMax, 1, yMin, yMax, 1, zMin, zMax
  )
  hist.SetStats(False)
  titleOffsets = (2.0, 2.0, 1.5)
  for index, axis in enumerate((hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis())):
    axis.CenterTitle(True)
    axis.SetTitleOffset(titleOffsets[index])
  hist.Draw()
  styleIndex = 0
  for dataSet, graph in parValueGraphs.items():
    graph.Draw("P ERR SAME")
    graph.SetName(dataSet)
    graph.SetTitle("")
    setCbFriendlyStyle(graph, styleIndex, skipBlack = False)
    styleIndex += 1
  legend = canv.BuildLegend()
  legend.SetFillStyle(0)
  legend.SetBorderSize(0)
  canv.SaveAs(f"{pdfDirName}/{hist.GetName()}.pdf")
  ROOT.gStyle.SetFrameFillColor(savedFrameFillColor)  # restore previous value


if __name__ == "__main__":
  printGitInfo()
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Plots BruFit results in given directory.")
  parser.add_argument("outputDirName", type = str, nargs = "?", default = "./BruFitOutput", help = "The path to the BruFit output directory; (default: '%(default)s')")
  args = parser.parse_args()
  dataSets    = ["Total", "Found", "Missing"]
  fitVariable = "MissingMassSquared_Measured"

  # plot fits and read parameter values from fit results
  parInfos:    Dict[str, List[ParInfo]]        = {}    # parInfos[<dataset>][<bin>]
  parNames:    Optional[Tuple[str, ...]]       = None
  binVarNames: Optional[List[Tuple[str, ...]]] = None  # binning variables for each binning
  for dataSet in dataSets:
    fitResultDirName  = f"{args.outputDirName}/{dataSet}"
    print(f"Plotting overall fit result for '{dataSet}' dataset")
    plotFitResult(BinInfo("", {}, {}, fitResultDirName), fitVariable)
    binVarNamesInDataSet: List[Tuple[str, ...]] = []
    for binningInfo in getBinningInfosFromDir(fitResultDirName):
      if binningInfo:
        binVarNamesInDataSet.append(binningInfo.varNames)
        if binningInfo.dirNames:
          print(f"Plotting fit results for binning variable(s) '{binningInfo.varNames}' for '{dataSet}' dataset")
          plotFitResults(binningInfo, fitVariable)
          if not dataSet in parInfos:
            parInfos[dataSet] = []
          parInfos[dataSet].extend(readParInfosForBinning(binningInfo))
          parNamesInBinning = parInfos[dataSet][-1].names  # readParInfosForBinning() ensures that parameter names are identical within a binning
          if parNames is not None:
            assert parNamesInBinning == parNames, f"The parameter set {parNamesInBinning} for dataset '{dataSet}' and binning '{binningInfo.varNames}' is different from the parameter set {parNames} of the previous dataset"
          else:
            parNames = parNamesInBinning
    if binVarNames is not None:
      assert binVarNamesInDataSet == binVarNames, f"The binning variables {binVarNamesInDataSet} for dataset '{dataSet}' are different from the binning variables {binVarNames} of the previous dataset"
    else:
      binVarNames = binVarNamesInDataSet

  # plot fit parameters as function of binning variable(s)
  setupPlotStyle()
  if parNames and binVarNames:
    for parName in parNames:
      for binningVars in binVarNames:
        if len(binningVars) == 1:
          plotParValue1D(parInfos, parName, binningVars[0],  args.outputDirName)
        elif len(binningVars) == 2:
          plotParValue2D(parInfos, parName, binningVars[:2], args.outputDirName)
