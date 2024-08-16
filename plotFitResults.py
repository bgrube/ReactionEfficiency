#!/usr/bin/env python3


from __future__ import annotations

import argparse
from collections import defaultdict
from collections.abc import (
  Generator,
  Iterator,
  Mapping,
  Sequence,
)
from dataclasses import dataclass
import functools
import itertools
import math
import os
import sys
from typing import Any

from uncertainties import UFloat, ufloat

import ROOT
if __name__ == "__main__":
  ROOT.PyConfig.DisableRootLogon = True  # do not change style of canvases loaded from fit result files

# import plotBeautifiers
from plotTools import (
  callMemberFunctionsWithArgs,
  drawZeroLine,
  getGraph1DFromValues,
  getGraph2DFromValues,
  printGitInfo,
  setCbFriendlyStyle,
  setupPlotStyle,
)


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


BINNING_VAR_PLOT_INFO: dict[str, dict[str, str]] = {
  "BeamEnergy"         : {"label" : "#it{E}_{beam}",                 "unit" : "GeV"},
  "MissingProtonP"     : {"label" : "#it{p}_{miss}^{kin. fit}",      "unit" : "GeV/#it{c}"},
  "MissingProtonTheta" : {"label" : "#it{#theta}_{miss}^{kin. fit}", "unit" : "deg"},
  "MissingProtonPhi"   : {"label" : "#it{#phi}_{miss}^{kin. fit}",   "unit" : "deg"},
}

def getAxisInfoForBinningVar(binningVar: str | tuple[str, str]) -> tuple[str, str, str]:
  """Returns name, label, and unit of the binning variable given by binningVar or binningVar[0]"""
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
  name:                str                             # bin name
  binDefs:             dict[str, tuple[float, float]]  # bin definitions { <binning var> : (<bin center>, <bin width>), ... }
  dirName:             str                             # directory name for bins
  nmbBootstrapSamples: int = 0                         # number of bootstrap samples

  def __len__(self) -> int:
    """Returns number of binning variables"""
    return len(self.binDefs)

  def __iter__(self) -> Generator[tuple[str, tuple[float, float]], None, None]:
    """Iterates over bin definitions"""
    return ((binName, (binDef[0], binDef[1])) for binName, binDef in self.binDefs.items())

  @property
  def centers(self) -> Generator[tuple[str, float], None, None]:
    """Returns generator for bin centers"""
    return ((binName, binDef[0]) for binName, binDef in self.binDefs.items())

  def center(
    self,
    varName: str,
  ) -> float:
    """Returns bin center for given binning variable"""
    if varName in self.binDefs:
      return self.binDefs[varName][0]
    else:
      raise KeyError(f"No bin center found for variable '{varName}'")

  @property
  def widths(self) -> Generator[tuple[str, float], None, None]:
    """Returns generator for bin widths"""
    return ((binName, binDef[1]) for binName, binDef in self.binDefs.items())

  def width(
    self,
    varName: str,
  ) -> float:
    """Returns bin width for given binning variable"""
    if varName in self.binDefs:
      return self.binDefs[varName][1]
    else:
      raise KeyError(f"No bin width found for variable '{varName}'")

  @property
  def varNames(self) -> tuple[str, ...]:
    """Returns names of kinematic variables used for binning"""
    return tuple(sorted(self.binDefs.keys()))

  @property
  def fitResultFileName(self) -> str:
    """Returns fit result file name for this bin"""
    return f"{self.dirName}/ResultsHSMinuit2.root"

  @property
  def bootstrapFileNames(self) -> Generator[tuple[int, str], None, None]:
    """Returns generator for file names with the bootstrap results for this bin"""
    return ((index, f"{self.dirName}/ResultsBoot{index}HSMinuit2.root") for index in range(self.nmbBootstrapSamples))

  def isSameBinAs(
    self,
    other: BinInfo,
  ) -> bool:
    """Returns whether 2 BinInfo objects represent the same kinematic bin"""
    return (self.name == other.name) and (self.binDefs == other.binDefs)


@dataclass
class BinningInfo:
  """Stores information about one particular kinematic binning"""
  infos:   list[BinInfo]  # infos for all bins of the kinematic binning
  dirName: str            # directory that contains all bins of the kinematic binning

  # make BinningInfo behave like a list of BinInfos
  def __len__(self) -> int:
    """Returns number of BinInfos"""
    return len(self.infos)

  def __getitem__(
    self,
    subscript: int
  ) -> BinInfo:
    """Returns BinInfo that correspond to given bin index"""
    return self.infos[subscript]

  def __iter__(self) -> Iterator[BinInfo]:
    """Iterates over BinInfos"""
    return iter(self.infos)

  @property
  def names(self) -> tuple[str, ...]:
    """Returns tuple with bin names"""
    return tuple(binInfo.name for binInfo in self.infos)

  @property
  def dirNames(self) -> tuple[str, ...]:
    """Returns tuple with directory names for all bins"""
    return tuple(binInfo.dirName for binInfo in self.infos)

  @property
  def fitResultFileNames(self) -> tuple[str, ...]:
    """Returns tuple with fit-result file names for all bin"""
    return tuple(binInfo.fitResultFileName for binInfo in self.infos)

  @property
  def varNames(self) -> tuple[str, ...]:
    """Returns names of kinematic variables used for binning"""
    varNames = None
    for binInfo in self.infos:
      varNamesInBin = binInfo.varNames
      if varNames is None:
        varNames = varNamesInBin
      else:
        assert varNamesInBin == varNames, f"Bins have inconsistent set of binning variables: bin '{binInfo.name}': {varNamesInBin} vs. {varNames} in previous bin"
    return varNames or tuple()

  def varBinCenters(
    self,
    varName: str,
  ) -> Generator[float, None, None]:
    """Returns generator for bin centers for given binning variable"""
    return (binInfo.center(varName) for binInfo in self.infos)

  def varBinWidths(
    self,
    varName: str,
  ) -> Generator[float, None, None]:
    """Returns generator for bin centers for given binning variable"""
    return (binInfo.width(varName) for binInfo in self.infos)

  def varRange(
    self,
    varName: str,
  ) -> tuple[float, float]:
    """Returns tuple with range of given binning variable"""
    # assumes non-overlapping bins
    varCenters  = [binInfo.center(varName) for binInfo in self.infos]
    minCenter   = min(varCenters)
    maxCenter   = max(varCenters)
    minIndex    = varCenters.index(minCenter)
    maxIndex    = varCenters.index(maxCenter)
    binWidthMin = self.infos[minIndex].width(varName)
    binWidthMax = self.infos[maxIndex].width(varName)
    return (minCenter - binWidthMin / 2.0, maxCenter + binWidthMax / 2.0)

  def varNmbBins(
    self,
    varName: str,
  ) -> int:
    """Returns number of bins for given binning variable"""
    # ensure equidistant bins
    widths = set(self.varBinWidths(varName))
    assert len(widths) == 1, f"Binning is not equidistant: '{varName}' bin widths = {widths}"
    width = next(iter(widths))  # get the only element
    varRange = self.varRange(varName)
    return int((varRange[1] - varRange[0]) / width)

  def varBinningInfo(
    self,
    varName: str,
  ) -> BinningInfo:
    """Returns BinningInfo object for 1-dimensional binning of given binning variable"""
    binInfos = [BinInfo(
      name    = binInfo.name,
      binDefs = {varName : (binInfo.center(varName), binInfo.width(varName))},
      dirName = binInfo.dirName,
    ) for binInfo in self.infos]
    return BinningInfo(binInfos, self.dirName)

  def isSameBinningAs(
    self,
    other: BinningInfo,
  ) -> bool:
    """Returns whether 2 BinningInfo objects represent the same kinematic binning"""
    if len(self.infos) != len(other.infos):
      return False
    for binInfo, otherBinInfo in zip(self.infos, other.infos):
      if not binInfo.isSameBinAs(otherBinInfo):
        raise ValueError("BinningInfo objects are not the same")
    return all(binInfo.isSameBinAs(otherBinInfo) for binInfo, otherBinInfo in zip(self.infos, other.infos))

  def __eq__(
    self,
    other: BinningInfo,
  ) -> bool:
    """Equality operator that compares binning information"""
    return self.isSameBinningAs(other)


def getBinningFromDir(fitResultDirName: str) -> BinningInfo | None:
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
  binInfos: list[BinInfo] = []
  for axisBinIndices in itertools.product(*axisBinIndexRanges):  # loop over all tuples of bin indices for the axes
    #TODO simplify code by moving this into dict comprehension below
    axisBinCenters = tuple(axes[axisIndex].GetBinCenter(axisBinIndex) for axisIndex, axisBinIndex in enumerate(axisBinIndices))
    axisBinWidths  = tuple(axes[axisIndex].GetBinWidth (axisBinIndex) for axisIndex, axisBinIndex in enumerate(axisBinIndices))
    binIndex = bins.FindBin(*axisBinCenters)  #!Note! the unpacking works only for up to 6 binning dimensions
    binName = binNames[binIndex]
    binInfos.append(BinInfo(
      name    = binName,
      binDefs = {binVarNames[axisIndex] : (axisBinCenter, axisBinWidth)
                 for axisIndex, (axisBinCenter, axisBinWidth) in enumerate(zip(axisBinCenters, axisBinWidths))},
      dirName = f"{fitResultDirName}/{str(binName)}",
    ))
  return BinningInfo(binInfos, fitResultDirName)


def getBinningInfosFromDir(fitResultDirName: str) -> list[BinningInfo | None]:
  """Reads binning infos from all 1st-level subdirectories that contain a 'DataBinsConfig.root' file"""
  binningInfos: list[BinningInfo | None] = []
  # find all subdirectories with binning files non-recursively
  subDirNames: list[str] = sorted([entry.path for entry in os.scandir(fitResultDirName) if entry.is_dir() and os.path.isfile(f"{entry.path}/DataBinsConfig.root")])
  # get binning info from all subdirectories
  for subDirName in subDirNames:
    print(f"Found binning info in directory '{subDirName}'")
    binningInfos.append(getBinningFromDir(subDirName))
  return binningInfos


@dataclass
class ParInfo:
  """Stores information about parameter values in a single kinematic bin"""
  binInfo:     BinInfo            # info for the kinematic bin
  fitVariable: str                # name of the fit variable
  values:      dict[str, UFloat]  # mapping of parameter names to values

  @property
  def names(self) -> tuple[str, ...]:
    """Returns parameter names"""
    return tuple(sorted(self.values.keys()))

  def readChi2From(
    self,
    fitResultFile: ROOT.TFile
  ) -> None:
    """Reads reduced chi^2 value from given fit-result file"""
    # BruFit (mis)uses TPad::SetTheta() to store the reduced chi^2
    # value in the canvas used to plot the fit result
    # see https://root.cern/doc/master/rf109__chi2residpull_8py.html
    # and https://root-forum.cern.ch/t/how-to-correctly-extract-the-chi2-ndf-p-value-of-a-roofitresult/45956
    canvName = (self.binInfo.name if self.binInfo.name and self.binInfo.name != "Overall" else "") + f"_{self.fitVariable}"
    canv = fitResultFile.Get(canvName)
    self.values["redChi2"] = ufloat(canv.GetTheta(), 0)


def acceptFitResult(fitResult: ROOT.RooFitResult) -> bool:
  """Returns whether the fit result should be accepted"""
  minimizerStatuses = [(i, fitResult.statusLabelHistory(i), fitResult.statusCodeHistory(i)) for i in range(fitResult.numStatusHistory())]
  if not all(status[2] == 0 for status in minimizerStatuses):
    print("Disregarding fit result because its status is " + ", ".join([f"[{status[0]}] {status[1]} = {status[2]}" for status in minimizerStatuses]))
    return False
  # if fitResult.covQual() != 3:  # this cut seems to be too aggressive, in particular when fudge parameters are freed
  #   print("Disregarding fit result because covariance matrix was not positive definite")
  #   return False
  if fitResult.minNll() == 0:
    print("Disregarding fit result because NLL at minimum = 0")
    return False
  if fitResult.edm() == 0:
    print("Disregarding fit result because EDM = 0")
    return False
  # remove fits with suspicious signal yields
  fitPars = fitResult.floatParsFinal()
  sigYieldIndex = fitPars.index("Yld_SigPdf")
  if sigYieldIndex < 0:
    print(f"Disregarding fit result because it does not contain the signal yield")
    return False
  sigYield = fitPars[sigYieldIndex]
  #TODO Use different thresholds for 'Found' and 'Missing' yields?
  if sigYield.getVal() < 100 or sigYield.getError() / sigYield.getVal() > 1:
    print(f"Disregarding fit result because the signal yield {sigYield.getVal()} +- {sigYield.getError()} is too small or its relative uncertainty too large")
    return False
  # remove fits with suspicious background yields
  bkgYieldIndex = fitPars.index("Yld_BkgPdf")
  if bkgYieldIndex < 0:
    print(f"Disregarding fit result because it does not contain the background yield")
    return False
  bkgYield = fitPars[bkgYieldIndex]
  if bkgYield.getVal() < 1:
    print(f"Disregarding fit result because the background yield {bkgYield.getVal()} is too small")
    return False
  #TODO cut in fit chi^2? problem is that this is stored in the canvas, not in the fit result; see readChi2From()
  return True


def readParInfoForBin(
  binInfo:           BinInfo,
  fitVariable:       str,  # name of the fit variable
  fitParNamesToRead: Mapping[str, str] | None = None,  # if { <new par name> : <par name>, ... } is set, only the given parameters are read, where `new par name` is the key used in the output
) -> ParInfo | None:
  """Reads parameter values from fit result for given kinematic bin; an optional mapping selects parameters to read and translates parameter names"""
  fitResultFileName = binInfo.fitResultFileName
  if not os.path.isfile(fitResultFileName):
    print(f"Cannot find file '{fitResultFileName}'. Skipping bin {binInfo}.")
    return ParInfo(binInfo, fitVariable, {})
  print(f"Reading fit result object 'MinuitResult' from file '{fitResultFileName}'")
  fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
  fitResult     = fitResultFile.Get("MinuitResult")
  # filter out unsuccessful fits
  if not acceptFitResult(fitResult):
    fitResultFile.Close()
    return None
  # print(f"{fitResult.numInvalidNLL()=}")
  # get fit parameters
  fitResult.Print()
  fitPars = fitResult.floatParsFinal()
  parValuesInBin: dict[str, UFloat] = {}
  if isinstance(fitParNamesToRead, Mapping):
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
  parInfo = ParInfo(binInfo, fitVariable, parValuesInBin)
  parInfo.readChi2From(fitResultFile)
  fitResultFile.Close()
  return parInfo


def readParInfosForBinning(
  binningInfo: BinningInfo,
  fitVariable: str,  # name of the fit variable
) -> list[ParInfo]:
  """Reads parameter values from fit results for given kinematic binning"""
  parInfos: list[ParInfo] = []
  parNames = None  # used to compare parameter names for current and previous bin
  for binInfo in binningInfo.infos:
    parInfo = readParInfoForBin(binInfo, fitVariable)
    if parInfo is not None and parInfo.values:
      # ensure that the parameter sets of the fit functions are the same in all kinematic bins, for which parameters could be read successfully
      parNamesInBin = parInfo.names
      if parNames is not None:
        assert parNamesInBin == parNames, f"The parameter set {parNamesInBin} of this bin is different from the parameter set {parNames} of the previous one"
      parNames = parNamesInBin
      print(f"Read parameter values for kinematic bin: {parInfo}")
      parInfos.append(parInfo)
  return parInfos


def plotFitResultForBin(
  binInfo:     BinInfo,
  fitVariable: str,  # name of the fit variable
  pdfDirName:  str | None = None,  # overrides default PDF output path (i.e. same dir as fit result file) if set
) -> None:
  """Plots fit result for given kinematic bin"""
  fitResultFileName = binInfo.fitResultFileName
  if not os.path.isfile(fitResultFileName):
    print(f"Cannot find file '{fitResultFileName}'. Skipping bin {binInfo}.")
    return
  optStatSave = ROOT.gStyle.GetOptStat()
  ROOT.gStyle.SetOptStat(False)
  print(f"Plotting fit result in file '{fitResultFileName}'")
  fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
  canvName = f"{binInfo.name or ''}_{fitVariable}"
  canv = fitResultFile.Get(canvName)
  dataFitPad = canv.GetListOfPrimitives().FindObject(f"{canvName}_1")
  # improve plot style
  curveStyles = {
    f"{binInfo.name or ''}TotalPDF_Norm[{fitVariable}]"              : {"SetLineColor" : ROOT.kRed   + 1, "SetLineWidth" : 2},
    f"{binInfo.name or ''}TotalPDF_Norm[{fitVariable}]_Comp[SigPdf]" : {"SetLineColor" : ROOT.kGreen + 2, "SetLineWidth" : 2},
    f"{binInfo.name or ''}TotalPDF_Norm[{fitVariable}]_Comp[BkgPdf]" : {"SetLineColor" : ROOT.kBlue  + 1, "SetLineWidth" : 2},
  }
  for curveName, style in curveStyles.items():
    curve = dataFitPad.GetListOfPrimitives().FindObject(curveName)
    callMemberFunctionsWithArgs(curve, style)
  hist = dataFitPad.GetListOfPrimitives().FindObject("h_DataEvents")
  global histClone
  histClone = hist.DrawClone("PZ")  # for some reason modifying draw options of hist has no effect; use copy instead
  histClone.SetLineColor(ROOT.kBlack)
  histClone.SetLineWidth(2)
  dataFitPad.GetListOfPrimitives().Remove(hist)
  # improve TPaveText with fit parameters
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
  del histClone
  ROOT.gStyle.SetOptStat(optStatSave)  # revert to previous setting


def plotFitResultsForBinning(
  binningInfo: BinningInfo,  # binning information
  fitVariable: str,  # name of the fit variable
) -> None:
  """Plots fit results for all kinematic bins"""
  for binInfo in binningInfo.infos:
    plotFitResultForBin(binInfo, fitVariable, pdfDirName = binningInfo.dirName)


def getParValuesForGraph1D(
  binVarName: str,  # name of the binning variable, i.e. x-axis
  parName:    str,  # name of value, i.e. y-axis
  parInfos:   Sequence[ParInfo],
) -> tuple[tuple[UFloat, UFloat], ...]:
  """Extracts information needed to plot parameter with given name as a function of the given bin variable from list of ParInfos"""
  graphValues: tuple[tuple[UFloat, UFloat], ...] = tuple(
    (ufloat(parInfo.binInfo.center(binVarName), parInfo.binInfo.width(binVarName) / 2.0), parInfo.values[parName])
    for parInfo in parInfos
    if (binVarName in parInfo.binInfo.varNames) and (parName in parInfo.names) and (len(parInfo.binInfo.varNames) == 1)
  )
  return graphValues


#TODO add sequence of beautifier functors as optional argument
def plotGraphs1D(
  graphOrGraphs:     ROOT.TGraph | Sequence[tuple[str, ROOT.TGraph]],
  binningInfo:       BinningInfo,  # 1D binning information
  yAxisTitle:        str,
  pdfDirName:        str,
  pdfFileBaseName:   str,
  pdfFileNamePrefix: str                        = "Proton_4pi_",
  pdfFileNameSuffix: str                        = "",
  graphTitle:        str | None                 = None,
  graphMinimum:      float | None               = 0.0,
  graphMaximum:      float | None               = None,
  skipBlack:         bool                       = True,
  drawLegend:        bool | None                = None,
  forceXRange:       tuple[float, float] | None = None,
  beautifiers:       Sequence[Any]              = [],  #TODO improve type hint by making a base class for beautifiers
) -> None:
  """Generic function that plots the given graph(s)"""
  assert len(binningInfo.varNames) == 1, f"Need 1-dimensional binning, but got {binningInfo}"
  binningVar = binningInfo.varNames[0]
  graphs: list[tuple[str, ROOT.TGraph]] = []
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
  #TODO make this function more generic by moving this into a beautifier
  _, binningVarLabel, binningVarUnit = getAxisInfoForBinningVar(binningVar)
  multiGraph.GetXaxis().SetTitle(f"{binningVarLabel} ({binningVarUnit})")
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
    multiGraph.GetXaxis().SetLimits(forceXRange[0], forceXRange[1])
  else:
    multiGraph.GetXaxis().SetLimits(*binningInfo.varRange(binningVar))
  for beautifier in beautifiers:
    beautifier.draw(multiGraph)
  if drawLegend == True or (drawLegend is None and legendLabels is not None):
    legend = ROOT.TLegend(0.3, 0.21)  # use automatic positioning
    for graph in multiGraph.GetListOfGraphs():
      label = graph.GetTitle()
      if not label:
        label = graph.GetName()
      legend.AddEntry(graph, label, "LPF")
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.Draw()
  # canv.Modified()
  # canv.Update()
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


def plotParValue1D(
  parInfos:          Mapping[str, Sequence[ParInfo]],
  parName:           str,  # name of parameter to plot
  binningInfo:       BinningInfo,  # 1D binning information
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str = "Proton_4pi_",
  pdfFileNameSuffix: str = "",
) -> None:
  """Overlays values of given parameter for given datasets for 1-dimensional binning with given binning variables"""
  assert len(binningInfo.varNames) == 1, f"Need 1-dimensional binning, but got {binningInfo}"
  binningVar = binningInfo.varNames[0]
  print(f"Plotting parameter '{parName}' as a function of binning variable '{binningVar}'")
  parValueGraphs: list[tuple[str, ROOT.TGraph]] = [(dataSet, getGraph1DFromValues(getParValuesForGraph1D(binningVar, parName, parInfosDataSet)))
                                                   for dataSet, parInfosDataSet in parInfos.items()]
  plotGraphs1D(
    graphOrGraphs     = parValueGraphs,
    binningInfo       = binningInfo,
    yAxisTitle        = parName,
    pdfDirName        = pdfDirName,
    pdfFileBaseName   = parName,
    pdfFileNamePrefix = pdfFileNamePrefix,
    pdfFileNameSuffix = pdfFileNameSuffix,
    graphTitle        = f"Fit parameter {parName})",
    skipBlack         = False,
  )


def getParValuesForGraph2D(
  binVarNames: Sequence[str],  # names of the binning variables, i.e. x-axis and y-axis
  parName:     str,  # name of value, i.e. z-axis
  parInfos:    Sequence[ParInfo],
) -> tuple[tuple[UFloat, UFloat, UFloat], ...]:
  """Extracts information needed to plot parameter with given name as a function of the given bin variables from list of ParInfos"""
  graphValues: tuple[tuple[UFloat, UFloat, UFloat], ...] = tuple(
    (
      ufloat(parInfo.binInfo.center(binVarNames[0]), parInfo.binInfo.width(binVarNames[0]) / 2.0),
      ufloat(parInfo.binInfo.center(binVarNames[1]), parInfo.binInfo.width(binVarNames[1]) / 2.0),
      parInfo.values[parName]
    )
    for parInfo in parInfos
    if (binVarNames[0] in parInfo.binInfo.varNames) and (binVarNames[1] in parInfo.binInfo.varNames) and (parName in parInfo.names) and (len(parInfo.binInfo.varNames) == 2)
  )
  return graphValues


def plotParValue2D(
  parInfos:          Mapping[str, Sequence[ParInfo]],
  parName:           str,  # name of parameter to plot
  binningInfo:       BinningInfo,  # 2D binning information
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str = "Proton_4pi_",
  pdfFileNameSuffix: str = "",
) -> None:
  """Overlays values of given parameter for given datasets for 2-dimensional binning with given binning variables"""
  assert len(binningInfo.varNames) == 2, f"Need 2-dimensional binning, but got {binningInfo}"
  binningVars = binningInfo.varNames[:2]
  savedFrameFillColor = ROOT.gStyle.GetFrameFillColor()
  ROOT.gStyle.SetFrameFillColor(0)  # switch back to default; otherwise graphs obstruct histogram frame
  print(f"Plotting parameter '{parName}' as a function of binning variables '{binningVars}'")
  parValueGraphs: dict[str, ROOT.TGraph2DErrors] = {}  # store graphs here to keep them in memory
  zMin = +math.inf
  zMax = -math.inf
  for dataSet in parInfos:
    graph = getGraph2DFromValues(getParValuesForGraph2D(binningVars, parName, parInfos[dataSet]))
    if graph is not None:
      parValueGraphs[dataSet] = graph
      zMin = min(zMin, graph.GetZmin())
      zMax = max(zMax, graph.GetZmax())
  if not parValueGraphs:  # nothing to plot
    return
  canv = ROOT.TCanvas()
  # construct dummy histogram that ensures same axes for all graphs
  binningVarLabels: list[str] = [""] * len(binningVars)
  binningVarUnits:  list[str] = [""] * len(binningVars)
  for index, binningVar in enumerate(binningVars):
    _, binningVarLabels[index], binningVarUnits[index] = getAxisInfoForBinningVar(binningVar)
  hist = ROOT.TH3F(
    f"{pdfFileNamePrefix}{parName}_{binningVars[0]}_{binningVars[1]}{pdfFileNameSuffix}",
    f";{binningVarLabels[0]} ({binningVarUnits[0]})"
    f";{binningVarLabels[1]} ({binningVarUnits[1]})"
    f";{parName}",
    1, *binningInfo.varRange(binningVars[0]),
    1, *binningInfo.varRange(binningVars[1]),
    1, zMin if zMin < 0 else 0, zMax,
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
  legend = ROOT.TLegend(0.3, 0.21)  # use automatic positioning
  for graph in parValueGraphs.values():
    legend.AddEntry(graph, graph.GetName(), "LPF")
  legend.SetFillStyle(0)
  legend.SetBorderSize(0)
  legend.Draw()
  canv.Update()  # needed to get rid of transient frame and axis objects
  canv.SaveAs(f"{pdfDirName}/{hist.GetName()}.pdf")
  ROOT.gStyle.SetFrameFillColor(savedFrameFillColor)  # restore previous value


def plotFitResults(
  dataSets:    Sequence[str] = ("Total", "Found", "Missing"),
  fitDirName:  str           = "./BruFitOutput",               # path to the BruFit output directory
  fitVariable: str           = "MissingMassSquared_Measured",  # name of the fit variable
) -> None:
  """Plots fit results and parameter values for given datasets in given directory"""
  parInfos:     defaultdict[str, list[ParInfo]] = defaultdict(list)  # parInfos[<dataset>][<bin>]
  parNames:     tuple[str, ...] | None          = None
  binningInfos: list[BinningInfo | None]        = []
  for dataSet in dataSets:
    fitResultDirName  = f"{fitDirName}/{dataSet}"
    print(f"Plotting overall fit result for '{dataSet}' dataset")
    plotFitResultForBin(BinInfo(name = "", binDefs = {}, dirName = fitResultDirName), fitVariable)
    binningInfosForDataSet: list[BinningInfo | None] = getBinningInfosFromDir(fitResultDirName)
    for binningInfo in binningInfosForDataSet:
      if binningInfo:
        # binningInfosForDataSet.append(binningInfo.varNames)
        if binningInfo.dirNames:
          print(f"Plotting fit results for binning variable(s) '{binningInfo.varNames}' for '{dataSet}' dataset")
          plotFitResultsForBinning(binningInfo, fitVariable)
          parInfos[dataSet].extend(readParInfosForBinning(binningInfo, fitVariable))
          parNamesInBinning = parInfos[dataSet][-1].names  # readParInfosForBinning() ensures that parameter names are identical within a binning
          if parNames is not None:
            assert parNamesInBinning == parNames, f"The parameter set {parNamesInBinning} for dataset '{dataSet}' and binning '{binningInfo.varNames}' is different from the parameter set {parNames} of the previous dataset"
          else:
            parNames = parNamesInBinning
    if binningInfos:
      assert len(binningInfos) == len(binningInfosForDataSet), f"The number of binnings {len(binningInfosForDataSet)} for dataset '{dataSet}' is different from the number of binnings {len(binningInfos)} of the previous dataset"
      for binningInfo, binningInfoDataSet in zip(binningInfos, binningInfosForDataSet):
        assert binningInfo == binningInfoDataSet, f"The binning {binningInfoDataSet} for dataset '{dataSet}' is different from the binning {binningInfo} of the previous dataset"
    else:
      binningInfos = binningInfosForDataSet

  # plot fit parameters as function of binning variable(s)
  if parNames and binningInfos:
    for parName in parNames:
      for binningInfo in binningInfos:
        if binningInfo:
          if len(binningInfo.varNames) == 1:
            plotParValue1D(parInfos, parName, binningInfo,  fitDirName)
          elif len(binningInfo.varNames) == 2:
            plotParValue2D(parInfos, parName, binningInfo, fitDirName)


if __name__ == "__main__":
  printGitInfo()
  ROOT.gROOT.SetBatch(True)
  setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Plots BruFit results in given directory.")
  parser.add_argument("fitDirName", type = str, nargs = "?", default = "./BruFitOutput", help = "The path to the BruFit output directory; (default: '%(default)s')")
  args = parser.parse_args()

  plotFitResults(fitDirName = args.fitDirName)
