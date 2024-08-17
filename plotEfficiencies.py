#!/usr/bin/env python3


from __future__ import annotations

import argparse
from collections import defaultdict
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

from plotFitResults import (
  BinInfo,
  BinningInfo,
  getAxisInfoForBinningVar,
  getBinningInfosFromDir,
  ParInfo,
  plotGraphs1D,
  readParInfoForBin,
)
from plotTools import (
  getGraph1DFromValues,
  getGraph2DFromValues,
  printGitInfo,
  redrawFrame,
  setupPlotStyle,
)


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


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


YIELD_PAR_NAMES: dict[str, str] = {
  "Signal"     : "Yld_SigPdf",
  "Background" : "Yld_BkgPdf",
}


def readYieldInfosForBinning(
  binningInfo:   BinningInfo,
  readIntegrals: bool = False,
  fitVariable:   str  = "",
) -> list[ParInfo]:
  """Reads yields (or histogram integrals if readIntegrals is set) from fit-result files for given binning"""
  yieldInfos: list[ParInfo] = []
  # read overall yields
  overallBinInfo = BinInfo(name = "Overall", binDefs = {}, dirName = binningInfo.dirName)  # special bin for overall fit results
  if os.path.isfile(overallBinInfo.fitResultFileName):
    yieldInfo = readDataIntegralFromFitFile(overallBinInfo, fitVariable) if readIntegrals \
      else readParInfoForBin(overallBinInfo, fitVariable, YIELD_PAR_NAMES)
    if yieldInfo is not None:
      print(f"Read overall yields: {yieldInfo}")
      yieldInfos.append(yieldInfo)  # first entry always contains overall yields
  # read yields for each bin
  for binInfo in binningInfo.infos:
    yieldInfo = readDataIntegralFromFitFile(binInfo, fitVariable) if readIntegrals \
      else readParInfoForBin(binInfo, fitVariable, YIELD_PAR_NAMES)
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
    (ufloat(effInfo.binInfo.center(binVarName), effInfo.binInfo.width(binVarName) / 2.0), effInfo.value)
    for effInfo in effInfos
    if (binVarName in effInfo.binInfo.varNames) and (len(effInfo.binInfo.varNames) == 1)]
  return graphValues


def plotEfficiencies1D(
  efficiencies:      Sequence[EffInfo],
  binningInfo:       BinningInfo,  # 1D binning information
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str = "Proton_4pi_",
  pdfFileNameSuffix: str = "",
) -> None:
  """Plots efficiency as a function of given binning variable for 1-dimensional binning"""
  assert len(binningInfo.varNames) == 1, f"Need 1-dimensional binning, but got {binningInfo}"
  binningVar = binningInfo.varNames[0]
  print(f"Plotting efficiency as a function of binning variable '{binningVar}'")
  plotGraphs1D(
    graphOrGraphs     = getGraph1DFromValues(getEffValuesForGraph1D(binningVar, efficiencies)),
    binningInfo       = binningInfo,
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
      ufloat(effInfo.binInfo.center(binVarNames[0]), effInfo.binInfo.width(binVarNames[0]) / 2.0),
      ufloat(effInfo.binInfo.center(binVarNames[1]), effInfo.binInfo.width(binVarNames[1]) / 2.0),
      effInfo.value
    )
    for effInfo in effInfos
    if (binVarNames[0] in effInfo.binInfo.varNames) and (binVarNames[1] in effInfo.binInfo.varNames) and (len(effInfo.binInfo.varNames) == 2)
  )
  return graphValues


def plotEfficiencies2DPerspective(
  efficiencies:      Sequence[EffInfo],
  binningInfo:       BinningInfo,  # 2D binning information
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str   = "Proton_4pi_",
  pdfFileNameSuffix: str   = "",
  markerSize:        float = 0.75,
) -> None:
  """Plots 3D view of efficiency as a function of given binning variables for 2-dimensional binning"""
  assert len(binningInfo.varNames) == 2, f"Need 2-dimensional binning, but got {binningInfo}"
  binningVars = binningInfo.varNames[:2]
  print(f"Plotting efficiency as a function of binning variables '{binningVars}'")
  efficiencyGraph = getGraph2DFromValues(getEffValuesForGraph2D(binningVars, efficiencies))
  if efficiencyGraph is None:  # nothing to plot
    return
  # ensure correct axis ranges using dummy histogram
  binningVarLabels: list[str] = [""] * len(binningVars)
  binningVarUnits:  list[str] = [""] * len(binningVars)
  for index, binningVar in enumerate(binningVars):
    _, binningVarLabels[index], binningVarUnits[index] = getAxisInfoForBinningVar(binningVar)
  hist = ROOT.TH2F(
    f"{pdfFileNamePrefix}mm2_eff_{binningVars[0]}_{binningVars[1]}{pdfFileNameSuffix}",
    f";{binningVarLabels[0]} ({binningVarUnits[0]})"
    f";{binningVarLabels[1]} ({binningVarUnits[1]})"
    ";Track-Finding Efficiency",
    1, *binningInfo.varRange(binningVars[0]),
    1, *binningInfo.varRange(binningVars[1]),
  )
  hist.SetStats(False)
  efficiencyGraph.SetHistogram(hist)
  efficiencyGraph.SetMinimum(0)
  efficiencyGraph.SetMaximum(1)
  canv = ROOT.TCanvas()
  efficiencyGraph.Draw("TRI2 P0 ERR")
  efficiencyGraph.SetTitle("")
  efficiencyGraph.SetMarkerStyle(ROOT.kFullCircle)
  efficiencyGraph.SetMarkerSize(markerSize)
  titleOffsets = (2.0, 2.0, 1.5)
  for index, axis in enumerate((efficiencyGraph.GetXaxis(), efficiencyGraph.GetYaxis(), efficiencyGraph.GetZaxis())):
    axis.CenterTitle(True)
    axis.SetTitleOffset(titleOffsets[index])
  # canv.Modified()
  # canv.Update()
  canv.SaveAs(f"{pdfDirName}/{hist.GetName()}.pdf")


def plotEfficiencies2DColzText(
  efficiencies:      Sequence[EffInfo],
  binningInfo:       BinningInfo,  # 2D binning information
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str   = "Proton_4pi_",
  pdfFileNameSuffix: str   = "",
) -> None:
  """Plots efficiency as a function of given binning variables for 2-dimensional binning using 'COLZ TEXT' option; works only for equidistant binning"""
  assert len(binningInfo.varNames) == 2, f"Need 2-dimensional binning, but got {binningInfo}"
  binningVars = tuple(reversed(binningInfo.varNames[:2]))  # swap p and theta axes
  print(f"Plotting efficiency as a function of binning variables '{binningVars}' assuming equidistant binning")
  # filter out relevant efficiencies
  effInfos = tuple(
    effInfo for effInfo in efficiencies
    if (binningVars[0] in effInfo.binInfo.varNames) and (binningVars[1] in effInfo.binInfo.varNames) and (len(effInfo.binInfo.varNames) == 2)
  )
  # TGraph2D always performs interpolation when drawn with COLZ -> construct TH2 with matching binning
  # create histogram
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}mm2_eff_{binningVars[0]}_{binningVars[1]}{pdfFileNameSuffix}", "")
  binningVarLabels: list[str] = [""] * len(binningVars)
  binningVarUnits:  list[str] = [""] * len(binningVars)
  for index, binningVar in enumerate(binningVars):
    _, binningVarLabels[index], binningVarUnits[index] = getAxisInfoForBinningVar(binningVar)
  efficiencyHist = ROOT.TH2D(
    f"h{canv.GetName()}",
    f";{binningVarLabels[0]} ({binningVarUnits[0]})"
    f";{binningVarLabels[1]} ({binningVarUnits[1]})",
    binningInfo.varNmbBins(binningVars[0]), *binningInfo.varRange(binningVars[0]),
    binningInfo.varNmbBins(binningVars[1]), *binningInfo.varRange(binningVars[1]),
  )
  # fill histogram
  for effInfo in effInfos:
    efficiencyHist.SetBinContent(efficiencyHist.FindBin(effInfo.binInfo.center(binningVars[0]), effInfo.binInfo.center(binningVars[1])),
                                 effInfo.value.nominal_value)
  # draw histogram
  efficiencyHist.SetMinimum(1e-9) # if set to exactly 0, zero entries are printed
  efficiencyHist.SetMaximum(1)
  ROOT.gStyle.SetPaintTextFormat("1.3f")
  efficiencyHist.Draw("COLZ TEXT")
  efficiencyHist.SetStats(False)
  redrawFrame(canv)
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}_ColzText.pdf")


def plotEfficiencies(
  dataSets:    Sequence[str] = ("Total", "Found", "Missing"),
  fitDirName:  str           = "./BruFitOutput",               # path to the BruFit output directory
  fitVariable: str           = "MissingMassSquared_Measured",  # name of the fit variable
  useMissing:  bool          = True,                           # True: 'Missing' dataset is used for efficiency calculation; False: 'Total' dataset is used instead
) -> None:
  """Plots efficiencies obtained from BruFit result in given directory"""
  for readIntegrals in (True, False):
    print("Calculating efficiencies from " + ("integrals of data histograms" if readIntegrals else "fit results"))
    yieldInfos:  defaultdict[str, list[ParInfo]] = defaultdict(list)    # [<dataset>][<bin>]
    binningInfos: list[BinningInfo | None]       = []
    for dataSet in dataSets:
      print("Reading " + ("integrals of data histograms" if readIntegrals else "yields") + f" for '{dataSet}' dataset")
      fitResultDirName = f"{fitDirName}/{dataSet}"
      yieldInfos[dataSet] = readYieldInfosForBinning(BinningInfo([], fitResultDirName), readIntegrals, fitVariable)  # overall yield values
      binningInfosForDataSet: list[BinningInfo | None] = getBinningInfosFromDir(fitResultDirName)
      for binningInfo in binningInfosForDataSet:
        if binningInfo:
          yieldInfos[dataSet].extend(readYieldInfosForBinning(binningInfo, readIntegrals, fitVariable))
      if binningInfos:
        assert len(binningInfos) == len(binningInfosForDataSet), f"The number of binnings {len(binningInfosForDataSet)} for dataset '{dataSet}' is different from the number of binnings {len(binningInfos)} of the previous dataset"
        for binningInfo, binningInfoDataSet in zip(binningInfos, binningInfosForDataSet):
          assert binningInfo == binningInfoDataSet, f"The binning {binningInfoDataSet} for dataset '{dataSet}' is different from the binning {binningInfo} of the previous dataset"
      else:
        binningInfos = binningInfosForDataSet

    effInfos = calculateEfficiencies(yieldInfos, useMissing)
    if effInfos:
      if effInfos[0].binInfo.name == "Overall":
        print("Overall efficiency from " + ("data-histogram integrals" if readIntegrals else "fits") + f" is {effInfos[0].value}")
      for binningInfo in binningInfos:
        if binningInfo:
          if len(binningInfo.varNames) == 1:
            plotEfficiencies1D(effInfos, binningInfo,  fitDirName, pdfFileNameSuffix = "_integral" if readIntegrals else "")
          elif len(binningInfo.varNames) == 2:
            plotEfficiencies2DPerspective(effInfos, binningInfo, fitDirName, pdfFileNameSuffix = "_integral" if readIntegrals else "")
            plotEfficiencies2DColzText   (effInfos, binningInfo, fitDirName, pdfFileNameSuffix = "_integral" if readIntegrals else "")


if __name__ == "__main__":
  printGitInfo()
  ROOT.gROOT.SetBatch(True)
  setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Plots efficiencies obtained from BruFit result in given directory.")
  parser.add_argument("--useTotal", action = "store_true", help = "If set, 'Total' distributions are used for efficiency calculation")
  parser.add_argument("fitDirName", type = str, nargs = "?", default = "./BruFitOutput", help = "The path to the BruFit output directory; (default: '%(default)s')")
  args = parser.parse_args()

  plotEfficiencies(fitDirName = args.fitDirName, useMissing = not args.useTotal)
