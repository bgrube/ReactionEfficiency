#!/usr/bin/env python3


from __future__ import annotations

import argparse
from collections.abc import (
  MutableMapping,
  Sequence,
)
from dataclasses import dataclass
import functools
import math
import os
import sys

from uncertainties import UFloat, ufloat

import ROOT
if __name__ == "__main__":
  ROOT.PyConfig.DisableRootLogon = True  # do not change style of canvases loaded from fit result files

import plotFitResults
from plotFitResults import (
  BinInfo,
  BinningInfo,
  getAxisInfoForBinningVar,
  ParInfo,
)
import plotTools


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


YIELD_PAR_NAMES: dict[str, str] = {
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
    return ParInfo(binInfo, fitVariable, {})
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
  return ParInfo(binInfo, fitVariable, {"Signal" : ufloat(integral, math.sqrt(integral))})


def readYieldInfosForBinning(
  binningInfo:   BinningInfo,
  readIntegrals: bool = False,
  fitVariable:   str  = "",
) -> list[ParInfo]:
  """Reads yields (or histogram integrals if readIntegrals is set) from fit-result files for given binning"""
  yieldInfos: list[ParInfo] = []
  # read overall yields
  overallBinInfo = BinInfo("Overall", {}, {}, binningInfo.dirName)  # special bin for overall fit results
  if os.path.isfile(overallBinInfo.fitResultFileName):
    yieldInfo = readDataIntegralFromFitFile(overallBinInfo, fitVariable) if readIntegrals \
      else plotFitResults.readParInfoForBin(overallBinInfo, fitVariable, YIELD_PAR_NAMES)
    if yieldInfo is not None:
      print(f"Read overall yields: {yieldInfo}")
      yieldInfos.append(yieldInfo)  # first entry always contains overall yields
  # read yields for each bin
  for binInfo in binningInfo.infos:
    yieldInfo = readDataIntegralFromFitFile(binInfo, fitVariable) if readIntegrals \
      else plotFitResults.readParInfoForBin(binInfo, fitVariable, YIELD_PAR_NAMES)
    if yieldInfo is not None:
      print(f"Read yields for kinematic bin: {yieldInfo}")
      yieldInfos.append(yieldInfo)
  return yieldInfos


@dataclass
class EffInfo:
  """Stores information about efficiency in a single kinematic bin"""
  binInfo: BinInfo  # info for the kinematic bin
  value:   UFloat   # efficiency value


def calculateEfficiencies(
  yieldInfos: MutableMapping[str, list[ParInfo]],  # [<dataset>][<bin>]
  useMissing: bool = True,  # True: 'Missing' dataset is used for efficiency calculation; False: 'Total' dataset is used instead
) -> list[EffInfo]:
  """Calculates efficiencies from yields"""
  # ensure yieldInfos contain the required datasets
  dataSets = ("Found", "Missing" if useMissing else "Total")
  for dataSet in dataSets:
    assert (dataSet in yieldInfos), f"'{dataSet}' dataset is missing"
  # ensure that yieldInfos for 'Found' and 'Missing'/'Total' have the same set of bins
  for setA, setB in (dataSets, reversed(dataSets)):
    yieldInfos[setA] = [
      yieldInfoFound for yieldInfoFound in yieldInfos[setA]
      if any(yieldInfoOther.binInfo.isSameBinAs(yieldInfoFound.binInfo) for yieldInfoOther in yieldInfos[setB])
    ]
  effInfos: list[EffInfo] = []
  for index, yieldInfoFound in enumerate(yieldInfos[dataSets[0]]):
    yieldInfoOther = yieldInfos[dataSets[1]][index]  # yield from 'Missing' or 'Total' dataset
    assert yieldInfoFound.binInfo.isSameBinAs(yieldInfoOther.binInfo), f"Bin infos for '{dataSets[0]}' and '{dataSets[1]}' are not identical: {yieldInfoFound.binInfo} vs. {yieldInfoOther.binInfo}"
    # calculate efficiency
    if "Signal" not in yieldInfoFound.names:
      print(f"No 'Signal' yield info for '{dataSets[0]}' case: {yieldInfoFound}")
      continue
    if "Signal" not in yieldInfoOther.names:
      print(f"No 'Signal' yield info for '{dataSets[1]}' case: {yieldInfoOther}")
      continue
    nmbFound = yieldInfoFound.values["Signal"]
    if useMissing:
      nmbTotal = nmbFound + yieldInfoOther.values["Signal"]  # 'Found' + 'Missing'
    else:
      nmbTotal = yieldInfoOther.values["Signal"]  # 'Total'
    effInfo = EffInfo(yieldInfoFound.binInfo, nmbFound / nmbTotal)
    print(f"Efficiency = {effInfo}")
    effInfos.append(effInfo)
  return effInfos


def getEffValuesForGraph1D(
  binVarName: str,  # name of the binning variable, i.e. x-axis
  effInfos:   Sequence[EffInfo],
) -> list[tuple[UFloat, UFloat]]:
  """Extracts information needed to plot efficiency as a function of the given bin variable from list of EffInfos"""
  graphValues: list[tuple[UFloat, UFloat]] = [
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
) -> tuple[tuple[UFloat, UFloat, UFloat], ...]:
  """Extracts information needed to plot efficiency as a function of the given bin variable from list of EffInfos"""
  graphValues: tuple[tuple[UFloat, UFloat, UFloat], ...] = tuple(
    (
      ufloat(effInfo.binInfo.centers[binVarNames[0]], effInfo.binInfo.widths[binVarNames[0]] / 2.0),
      ufloat(effInfo.binInfo.centers[binVarNames[1]], effInfo.binInfo.widths[binVarNames[1]] / 2.0),
      effInfo.value
    )
    for effInfo in effInfos
    if (binVarNames[0] in effInfo.binInfo.varNames) and (binVarNames[1] in effInfo.binInfo.varNames) and (len(effInfo.binInfo.varNames) == 2)
  )
  return graphValues


def plotEfficiencies2DPerspective(
  efficiencies:      Sequence[EffInfo],
  binningVars:       Sequence[str],  # names of binning variables to plot
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str   = "Proton_4pi_",
  pdfFileNameSuffix: str   = "",
  markerSize:        float = 0.75,
) -> None:
  """Plots 3D view of efficiency as a function of given binning variables for 2-dimensional binning"""
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
  binningVarLabels: list[str] = [""] * len(binningVars)
  binningVarUnits:  list[str] = [""] * len(binningVars)
  for index, binningVar in enumerate(binningVars):
    _, binningVarLabels[index], binningVarUnits[index] = getAxisInfoForBinningVar(binningVar)
  axisTitles = (
    f"{binningVarLabels[0]} ({binningVarUnits[0]})",
    f"{binningVarLabels[1]} ({binningVarUnits[1]})",
    f"Track-Finding Efficiency",
  )
  titleOffsets = (2.0, 2.0, 1.5)
  for index, axis in enumerate((efficiencyGraph.GetXaxis(), efficiencyGraph.GetYaxis(), efficiencyGraph.GetZaxis())):
    axis.SetTitle(axisTitles[index])
    axis.CenterTitle(True)
    axis.SetTitleOffset(titleOffsets[index])
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


def plotEfficiencies2DColzText(
  efficiencies:      Sequence[EffInfo],
  binningVarsIn:     Sequence[str],  # names of binning variables to plot
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str   = "Proton_4pi_",
  pdfFileNameSuffix: str   = "",
) -> None:
  """Plots efficiency as a function of given binning variables for 2-dimensional binning using 'COLZ TEXT' option; works only for equidistant binning"""
  binningVars = tuple(reversed(binningVarsIn))  # swap p and theta axes
  print(f"Plotting efficiency as a function of binning variables '{binningVars}' assuming equidistant binning")
  # filter out relevant efficiencies
  effInfos = tuple(effInfo for effInfo in efficiencies
                   if (binningVars[0] in effInfo.binInfo.varNames) and (binningVars[1] in effInfo.binInfo.varNames) and (len(effInfo.binInfo.varNames) == 2))
  # TGraph2D always performs interpolation when drawn with COLZ -> construct TH2 with matching binning
  # determine equidistant binning and create histogram
  xWidths = set(effInfo.binInfo.widths[binningVars[0]] for effInfo in effInfos)
  yWidths = set(effInfo.binInfo.widths[binningVars[1]] for effInfo in effInfos)
  assert (len(xWidths) == 1) and (len(yWidths) == 1), f"Binning is not equidistant: x bin widths = {xWidths}; y bin widths = {yWidths}"
  xWidth = tuple(xWidths)[0]
  yWidth = tuple(yWidths)[0]
  xCenters = sorted(set(effInfo.binInfo.centers[binningVars[0]] for effInfo in effInfos))
  yCenters = sorted(set(effInfo.binInfo.centers[binningVars[1]] for effInfo in effInfos))
  xRange = (xCenters[0] - xWidth / 2.0, xCenters[-1] + xWidth / 2.0)
  yRange = (yCenters[0] - yWidth / 2.0, yCenters[-1] + yWidth / 2.0)
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}mm2_eff_{binningVars[0]}_{binningVars[1]}{pdfFileNameSuffix}", "")
  binningVarLabels: list[str] = [""] * len(binningVars)
  binningVarUnits:  list[str] = [""] * len(binningVars)
  for index, binningVar in enumerate(binningVars):
    _, binningVarLabels[index], binningVarUnits[index] = getAxisInfoForBinningVar(binningVar)
  efficiencyHist = ROOT.TH2D(
    f"h{canv.GetName()}", f";{binningVarLabels[0]} ({binningVarUnits[0]})"
                          f";{binningVarLabels[1]} ({binningVarUnits[1]})",
    int((xRange[1] - xRange[0]) / xWidth), *xRange, int((yRange[1] - yRange[0]) / yWidth), *yRange)
  # fill histogram
  for effInfo in effInfos:
    efficiencyHist.SetBinContent(efficiencyHist.FindBin(effInfo.binInfo.centers[binningVars[0]], effInfo.binInfo.centers[binningVars[1]]),
                                 effInfo.value.nominal_value)
  # draw histogram
  efficiencyHist.SetMinimum(1e-9) # if set to exactly 0, zero entries are printed
  efficiencyHist.SetMaximum(1)
  ROOT.gStyle.SetPaintTextFormat("1.3f")
  efficiencyHist.Draw("COLZ TEXT")
  efficiencyHist.SetStats(False)
  plotTools.redrawFrame(canv)
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}_ColzText.pdf")


if __name__ == "__main__":
  plotTools.printGitInfo()
  ROOT.gROOT.SetBatch(True)
  plotTools.setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Plots efficiencies obtained from BruFit result in given directory.")
  parser.add_argument("--useTotal",    action = "store_true", help = "If set, 'Total' distributions are used for efficiency calculation")
  parser.add_argument("outputDirName", type = str, nargs = "?", default = "./BruFitOutput", help = "The path to the BruFit output directory; (default: '%(default)s')")
  args = parser.parse_args()
  dataSets    = ["Total", "Found", "Missing"]
  fitVariable = "MissingMassSquared_Measured"

  for readIntegrals in (True, False):
    print("Calculating efficiencies from " + ("integrals of data histograms" if readIntegrals else "fit results"))
    yieldInfos:  dict[str, list[ParInfo]]     = {}    # [<dataset>][<bin>]
    binVarNames: list[tuple[str, ...]] | None = None  # binning variables for each binning
    for dataSet in dataSets:
      print("Reading " + ("integrals of data histograms" if readIntegrals else "yields") + f" for '{dataSet}' dataset")
      fitResultDirName = f"{args.outputDirName}/{dataSet}"
      yieldInfos[dataSet] = readYieldInfosForBinning(BinningInfo([], fitResultDirName), readIntegrals, fitVariable)  # overall yield values
      binVarNamesInDataSet: list[tuple[str, ...]] = []
      for binningInfo in plotFitResults.getBinningInfosFromDir(fitResultDirName):
        if binningInfo:
          binVarNamesInDataSet.append(binningInfo.varNames)
          yieldInfos[dataSet].extend(readYieldInfosForBinning(binningInfo, readIntegrals, fitVariable))
      if binVarNames is None:
        binVarNames = binVarNamesInDataSet
      else:
        assert binVarNamesInDataSet == binVarNames, f"The binning variables {binVarNamesInDataSet} for dataset '{dataSet}' are different from the binning variables {binVarNames} of the previous dataset"
    effInfos = calculateEfficiencies(yieldInfos, useMissing = not args.useTotal)
    if effInfos:
      if effInfos[0].binInfo.name == "Overall":
        print("Overall efficiency from " + ("data-histogram integrals" if readIntegrals else "fits") + f" is {effInfos[0].value}")
      if binVarNames:
        for binningVars in binVarNames:
          if len(binningVars) == 1:
            plotEfficiencies1D(effInfos, binningVars[0],  args.outputDirName, pdfFileNameSuffix = "_integral" if readIntegrals else "")
          elif len(binningVars) == 2:
            plotEfficiencies2DPerspective(effInfos, binningVars[:2], args.outputDirName, pdfFileNameSuffix = "_integral" if readIntegrals else "")
            plotEfficiencies2DColzText   (effInfos, binningVars[:2], args.outputDirName, pdfFileNameSuffix = "_integral" if readIntegrals else "")
