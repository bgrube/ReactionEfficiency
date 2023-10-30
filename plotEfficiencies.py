#!/usr/bin/env python3


import argparse
from dataclasses import dataclass
import functools
import math
import numpy as np
import os
import sys
from typing import (
  Dict,
  List,
  Mapping,
  Optional,
  Sequence,
  Tuple,
)

from uncertainties import UFloat, ufloat

import ROOT

import makePlots
import plotFitResults
from plotFitResults import (
  BinInfo,
  BinningInfo,
  ParInfo,
  BINNING_VAR_PLOT_INFO,
)
import plotTools


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


YIELD_PAR_NAMES: Dict[str, str] = {
  "Signal"     : "Yld_SigPdf",
  "Background" : "Yld_BkgPdf",
}


def readDataIntegralFromFitFile(
  binInfo:     BinInfo,
  fitVariable: str,
) -> ParInfo:
  """Reads data histogram for given kinematic bin and returns its integral assuming Poissonian uncertainty"""
  assert fitVariable, f"Cannot use value '{fitVariable}' as fit variable name"
  fitResultFileName = binInfo.fitResultFileName
  if not os.path.isfile(fitResultFileName):
    print(f"Cannot find file '{fitResultFileName}'. Skipping bin {binInfo}.")
    return ParInfo(binInfo, {})
  print(f"Reading fitted data histogram from file '{fitResultFileName}'")
  fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
  canvName = f"{'' if binInfo.name == 'Overall' else binInfo.name}_{fitVariable}"
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
  return ParInfo(binInfo, {"Signal" : ufloat(integral, math.sqrt(integral))})


def readYieldInfosForBinning(
  binningInfo:   BinningInfo,
  readIntegrals: bool = False,
  fitVariable:   str  = "",
) -> List[ParInfo]:
  """Reads yields (or histogram integrals if readIntegrals is set) from fit-result files for given binning"""
  yieldInfos = []
  # read overall yields
  overallBinInfo = BinInfo("Overall", {}, {}, binningInfo.dirName)  # special bin for overall fit results
  if os.path.isfile(overallBinInfo.fitResultFileName):
    yieldInfo = readDataIntegralFromFitFile(overallBinInfo, fitVariable) if readIntegrals \
      else plotFitResults.readParInfoForBin(overallBinInfo, YIELD_PAR_NAMES)
    print(f"Read overall yields: {yieldInfo}")
    yieldInfos.append(yieldInfo)  # first entry always contains overall yields
  # read yields for each bin
  for binInfo in binningInfo.infos:
    yieldInfo = readDataIntegralFromFitFile(binInfo, fitVariable) if readIntegrals \
      else plotFitResults.readParInfoForBin(binInfo, YIELD_PAR_NAMES)
    print(f"Read yields for kinematic bin: {yieldInfo}")
    yieldInfos.append(yieldInfo)
  return yieldInfos


@dataclass
class EffInfo:
  """Stores information about efficiency in a single kinematic bin"""
  binInfo: BinInfo  # info for the kinematic bin
  value:   UFloat   # efficiency value


def calculateEfficiencies(yieldInfos: Mapping[str, Sequence[ParInfo]]) -> List[EffInfo]:
  """Calculates efficiencies from yields"""
  assert ("Found" in yieldInfos) and ("Missing" in yieldInfos), "Either 'Found', 'Missing' or both datasets are missing"
  assert len(yieldInfos["Found"]) == len(yieldInfos["Missing"]), f"'Found' and 'Missing' datasets have different number of kinematic bins: {len(yieldInfos['Found'])} vs. {len(yieldInfos['Missing'])}"
  effInfos = []
  for index, yieldInfoFound in enumerate(yieldInfos["Found"]):
    yieldInfoMissing = yieldInfos["Missing"][index]
    assert yieldInfoFound.binInfo.isSameBinAs(yieldInfoMissing.binInfo), f"Bin infos for 'Found' and 'Missing' are not identical: {yieldInfoFound.binInfo} vs. {yieldInfoMissing.binInfo}"
    # calculate efficiency
    if "Signal" not in yieldInfoFound.names:
      print(f"No 'Signal' yield info for 'Found' case: {yieldInfoFound}")
      continue
    if "Signal" not in yieldInfoMissing.names:
      print(f"No 'Signal' yield info for 'Missing' case: {yieldInfoMissing}")
      continue
    nmbFound   = yieldInfoFound.values  ["Signal"]
    nmbMissing = yieldInfoMissing.values["Signal"]
    effInfo = EffInfo(yieldInfoFound.binInfo, nmbFound / (nmbFound + nmbMissing))
    print(f"Efficiency = {effInfo}")
    effInfos.append(effInfo)
  return effInfos


def getEffValuesForGraph1D(
  binVarName: str,  # name of the binning variable, i.e. x-axis
  effInfos:   Sequence[EffInfo],
) -> List[Tuple[UFloat, UFloat]]:
  """Extracts information needed to plot efficiency as a function of the given bin variable from list of EffInfos"""
  graphValues: List[Tuple[UFloat, UFloat]] = [
    (ufloat(effInfo.binInfo.centers[binVarName], effInfo.binInfo.widths[binVarName] / 2.0), effInfo.value)
    for effInfo in effInfos
    if (binVarName in effInfo.binInfo.varNames) and (len(effInfo.binInfo.varNames) == 1)]
  return graphValues


def plotEfficiencies1D(
  efficiencies:      Sequence[EffInfo],
  binningVar:        str,  # name of binning variable to plot
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str = "Proton_4pi_",
  pdfFileNameSuffix: str = "",
) -> None:
  """Plots efficiency as a function of given binning variable for 1-dimensional binning"""
  print(f"Plotting efficiency as a function of binning variable '{binningVar}'")
  plotFitResults.plotGraphs1D(
    graphOrGraphs     = plotTools.getGraph1DFromValues(getEffValuesForGraph1D(binningVar, efficiencies)),
    binningVar        = binningVar,
    yAxisTitle        = f"Track-Finding Efficiency",
    pdfDirName        = pdfDirName,
    pdfFileBaseName   = "mm2_eff",
    pdfFileNamePrefix = pdfFileNamePrefix,
    pdfFileNameSuffix = pdfFileNameSuffix,
    graphMinimum      = 0.0,
    graphMaximum      = 1.0,
    skipBlack         = False,
  )


def getEffValuesForGraph2D(
  binVarNames: Sequence[str],  # names of the binning variables, i.e. x-axis and y-axis
  effInfos:    Sequence[EffInfo],
) -> Tuple[Tuple[UFloat, UFloat, UFloat], ...]:
  """Extracts information needed to plot efficiency as a function of the given bin variable from list of EffInfos"""
  graphValues: Tuple[Tuple[UFloat, UFloat, UFloat], ...] = tuple(
    (
      ufloat(effInfo.binInfo.centers[binVarNames[0]], effInfo.binInfo.widths[binVarNames[0]] / 2.0),
      ufloat(effInfo.binInfo.centers[binVarNames[1]], effInfo.binInfo.widths[binVarNames[1]] / 2.0),
      effInfo.value
    )
    for effInfo in effInfos
    if (binVarNames[0] in effInfo.binInfo.varNames) and (binVarNames[1] in effInfo.binInfo.varNames) and (len(effInfo.binInfo.varNames) == 2)
  )
  return graphValues


def plotEfficiencies2D(
  efficiencies:      Sequence[EffInfo],
  binningVars:       Sequence[str],  # names of binning variables to plot
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str   = "Proton_4pi_",
  pdfFileNameSuffix: str   = "",
  markerSize:        float = 0.75,
) -> None:
  """Plots efficiency as a function of given binning variables for 2-dimensional binning"""
  print(f"Plotting efficiency as a function of binning variables '{binningVars}'")
  efficiencyGraph = plotTools.getGraph2DFromValues(getEffValuesForGraph2D(binningVars, efficiencies))
  if efficiencyGraph is None:  # nothing to plot
    return
  efficiencyGraph.SetMinimum(0)
  efficiencyGraph.SetMaximum(1)
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}mm2_eff_{binningVars[0]}_{binningVars[1]}{pdfFileNameSuffix}", "")
  efficiencyGraph.Draw("TRI2 P0 ERR")
  efficiencyGraph.SetTitle("")
  efficiencyGraph.SetMarkerStyle(ROOT.kFullCircle)
  efficiencyGraph.SetMarkerSize(markerSize)
  assert binningVars[0] in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVars[0]}'"
  assert binningVars[1] in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVars[1]}'"
  axisTitles = (
    f"{BINNING_VAR_PLOT_INFO[binningVars[0]]['label']} ({BINNING_VAR_PLOT_INFO[binningVars[0]]['unit']})",
    f"{BINNING_VAR_PLOT_INFO[binningVars[1]]['label']} ({BINNING_VAR_PLOT_INFO[binningVars[1]]['unit']})",
    f"Track-Finding Efficiency",
  )
  titleOffsets = (2.0, 2.0, 1.5)
  for index, axis in enumerate((efficiencyGraph.GetXaxis(), efficiencyGraph.GetYaxis(), efficiencyGraph.GetZaxis())):
    axis.SetTitle(axisTitles[index])
    axis.CenterTitle(True)
    axis.SetTitleOffset(titleOffsets[index])
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


if __name__ == "__main__":
  makePlots.printGitInfo()
  ROOT.gROOT.SetBatch(True)
  plotTools.setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Plots efficiencies obtained from BruFit result in given directory.")
  parser.add_argument("outputDirName", type = str, nargs = "?", default = "./BruFitOutput", help = "The path to the BruFit output directory; (default: '%(default)s')")
  args = parser.parse_args()
  dataSets    = ["Total", "Found", "Missing"]
  fitVariable = "MissingMassSquared_Measured"

  for readIntegrals in (True, False):
    print("Calculating efficiencies from " + ("integrals of data histograms" if readIntegrals else "fit results"))
    yieldInfos:  Dict[str, List[ParInfo]]        = {}    # yieldInfos[<dataset>][<bin>]
    binVarNames: Optional[List[Tuple[str, ...]]] = None  # binning variables for each binning
    for dataSet in dataSets:
      print("Reading " + ("integrals of data histograms" if readIntegrals else "yields") + f" for '{dataSet}' dataset")
      fitResultDirName = f"{args.outputDirName}/{dataSet}"
      yieldInfos[dataSet] = readYieldInfosForBinning(BinningInfo([], fitResultDirName), readIntegrals, fitVariable)  # overall yield values
      binVarNamesInDataSet = []
      for binningInfo in plotFitResults.getBinningInfosFromDir(fitResultDirName):
        if binningInfo:
          binVarNamesInDataSet.append(binningInfo.varNames)
          yieldInfos[dataSet].extend(readYieldInfosForBinning(binningInfo, readIntegrals, fitVariable))
      if binVarNames is None:
        binVarNames = binVarNamesInDataSet
      else:
        assert binVarNamesInDataSet == binVarNames, f"The binning variables {binVarNamesInDataSet} for dataset '{dataSet}' are different from the binning variables {binVarNames} of the previous dataset"
    effInfos = calculateEfficiencies(yieldInfos)
    if effInfos:
      if effInfos[0].binInfo.name == "Overall":
        print("Overall efficiency from " + ("data-histogram integrals" if readIntegrals else "fits") + f" is {effInfos[0].value}")
      if binVarNames:
        for binningVars in binVarNames:
          if len(binningVars) == 1:
            plotEfficiencies1D(effInfos, binningVars[0],  args.outputDirName, pdfFileNameSuffix = "_integral" if readIntegrals else "")
          elif len(binningVars) == 2:
            plotEfficiencies2D(effInfos, binningVars[:2], args.outputDirName, pdfFileNameSuffix = "_integral" if readIntegrals else "")
