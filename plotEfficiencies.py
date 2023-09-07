#!/usr/bin/env python3


import argparse
from dataclasses import dataclass  # builtin in Python 3.7+
import functools
import math
import numpy as np
import os
import sys
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from uncertainties import UFloat, ufloat

import ROOT

import makePlots
import plotFitResults
from plotFitResults import BinInfo, BinningInfo, ParInfo, BINNING_VAR_PLOT_INFO


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


YIELD_PAR_NAMES: Dict[str, str] = {
  "Signal"     : "Yld_SigPdf",
  "Background" : "Yld_BkgPdf",
}


def readDataIntegralFromFitFile(
  binInfo:     BinInfo,
  fitVariable: str,
) -> Optional[ParInfo]:
  """Reads data histogram for given kinematic bin and returns its integral assuming Poissonian uncertainty"""
  assert fitVariable, f"Cannot use value '{fitVariable}' as fit variable name"
  fitResultFileName = binInfo.fitResultFileName
  if not os.path.isfile(fitResultFileName):
    print(f"Cannot find file '{fitResultFileName}'. Skipping bin {binInfo}.")
    return None
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
    if yieldInfo is not None:
      print(f"Read overall yields: {yieldInfo}")
      yieldInfos.append(yieldInfo)  # first entry always contains overall yields
  # read yields for each bin
  for binInfo in binningInfo.infos:
    yieldInfo = readDataIntegralFromFitFile(binInfo, fitVariable) if readIntegrals \
      else plotFitResults.readParInfoForBin(binInfo, YIELD_PAR_NAMES)
    if yieldInfo is not None:
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
    effInfo = EffInfo(yieldInfoFound.binInfo,
      yieldInfoFound.values["Signal"] / (yieldInfoFound.values["Signal"] + yieldInfoMissing.values["Signal"]))
    print(f"Efficiency = {effInfo}")
    effInfos.append(effInfo)
  return effInfos


def getEffValuesForGraph1D(
  binVarName: str,  # name of the binning variable, i.e. x-axis
  effInfos:   Sequence[EffInfo],
) -> List[Tuple[float, UFloat]]:
  """Extracts information needed to plot efficiency as a function of the given bin variable from list of EffInfos"""
  graphValues: List[Tuple[float, UFloat]] = [(effInfo.binInfo.centers[binVarName], effInfo.value)
    for effInfo in effInfos if binVarName in effInfo.binInfo.varNames]
  return graphValues


def plotEfficiencies1D(
  efficiencies:      Sequence[EffInfo],
  binningVar:        str,  # name of binning variable to plot
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNameSuffix: str   = "",
  particle:          str   = "Proton",
  channel:           str   = "4pi",
  markerSize:        float = 0.75,
) -> None:
  """Plots efficiency as a function of given binning variable for 1-dimensional binning"""
  print(f"Plotting efficiency as a function of binning variable '{binningVar}'")
  efficiencyGraph = plotFitResults.getParValueGraph1D(getEffValuesForGraph1D(binningVar, efficiencies))
  # efficiencyGraph.SetTitle(f"{particle} Track-Finding Efficiency ({channel})")
  efficiencyGraph.SetMarkerStyle(ROOT.kFullCircle)
  efficiencyGraph.SetMarkerSize(markerSize)
  assert binningVar in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVar}'"
  efficiencyGraph.GetXaxis().SetTitle(f"{BINNING_VAR_PLOT_INFO[binningVar]['label']} ({BINNING_VAR_PLOT_INFO[binningVar]['unit']})")
  efficiencyGraph.GetYaxis().SetTitle(f"{particle} Track-Finding Efficiency")
  efficiencyGraph.SetMinimum(0)
  efficiencyGraph.SetMaximum(1)
  canv = ROOT.TCanvas(f"{particle}_{channel}_mm2_eff_{binningVar}{pdfFileNameSuffix}", "")
  efficiencyGraph.Draw("APZ")
  # #TODO? # indicate value from fit of overall distributions
  # line = ROOT.TLine()
  # line.SetLineStyle(ROOT.kDashed)
  # line.DrawLine(efficienciesKinBinsGraph.GetXaxis().GetXmin(), overallEff.nominal_value, efficienciesKinBinsGraph.GetXaxis().GetXmax(), overallEff.nominal_value)
  # # indicate weighted average of efficiencies in kinematic bins
  # meanEff = np.average(yVals, weights = [1 / (yErr**2) for yErr in yErrs])
  # line.SetLineColor(ROOT.kRed + 1)
  # line.DrawLine(efficiencyGraph.GetXaxis().GetXmin(), meanEff, efficiencyGraph.GetXaxis().GetXmax(), meanEff)
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


def getEffValuesForGraph2D(
  binVarNames: Sequence[str],  # names of the binning variables, i.e. x-axis and y-axis
  effInfos:    Sequence[EffInfo],
) -> Tuple[Tuple[float, float, UFloat], ...]:
  """Extracts information needed to plot efficiency as a function of the given bin variable from list of EffInfos"""
  graphValues: Tuple[Tuple[float, float, UFloat], ...] = tuple(
    (effInfo.binInfo.centers[binVarNames[0]], effInfo.binInfo.centers[binVarNames[1]], effInfo.value)
    for effInfo in effInfos if (binVarNames[0] in effInfo.binInfo.varNames) and (binVarNames[1] in effInfo.binInfo.varNames)
  )
  return graphValues


def plotEfficiencies2D(
  efficiencies:      Sequence[EffInfo],
  binningVars:       Sequence[str],  # names of binning variables to plot
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNameSuffix: str   = "",
  particle:          str   = "Proton",
  channel:           str   = "4pi",
  markerSize:        float = 0.75,
) -> None:
  """Plots efficiency as a function of given binning variables for 2-dimensional binning"""
  print(f"Plotting efficiency as a function of binning variables '{binningVars}'")
  efficiencyGraph = plotFitResults.getParValueGraph2D(getEffValuesForGraph2D(binningVars, efficiencies))
  if efficiencyGraph is None:  # nothing to plot
    return
  efficiencyGraph.SetMinimum(0)
  efficiencyGraph.SetMaximum(1)
  canv = ROOT.TCanvas(f"{particle}_{channel}_mm2_eff_{binningVars[0]}_{binningVars[1]}_{pdfFileNameSuffix}", "")
  efficiencyGraph.Draw("TRI2 P0")
  efficiencyGraph.SetTitle("")
  efficiencyGraph.SetMarkerStyle(ROOT.kFullCircle)
  efficiencyGraph.SetMarkerSize(markerSize)
  assert binningVars[0] in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVars[0]}'"
  assert binningVars[1] in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVars[1]}'"
  axisTitles = (
    f"{BINNING_VAR_PLOT_INFO[binningVars[0]]['label']} ({BINNING_VAR_PLOT_INFO[binningVars[0]]['unit']})",
    f"{BINNING_VAR_PLOT_INFO[binningVars[1]]['label']} ({BINNING_VAR_PLOT_INFO[binningVars[1]]['unit']})",
    f"{particle} Track-Finding Efficiency",
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
  makePlots.setupPlotStyle()
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
            plotEfficiencies1D(effInfos, binningVars[0],  args.outputDirName, "_integral" if readIntegrals else "")
          elif len(binningVars) == 2:
            plotEfficiencies2D(effInfos, binningVars[:2], args.outputDirName, "_integral" if readIntegrals else "")
