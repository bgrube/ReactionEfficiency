#!/usr/bin/env python3


import functools
import os
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
import overlayEfficiencies
import plotEfficiencies
from plotEfficiencies import EffInfo, BinInfo
import plotFitResults
from plotFitResults import BINNING_VAR_PLOT_INFO
import plotTools


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


def overlayEfficiencyRatios1D(
  effInfos:          Mapping[str, Mapping[Tuple[str, str], Sequence[EffInfo]]],
  binningVar:        str,
  pdfDirName:        str,
  pdfFileNameSuffix: str = "",
  particle:          str = "Proton",
  channel:           str = "4pi",
  graphTitle:        Optional[str] = None,
) -> None:
  """Plots efficiency ratios as a function of `binningVar` for all given fits with 1D binning"""
  print(f"Plotting efficiency ratio for binning variable '{binningVar}'")
  ratioGraphs: List[Tuple[str, ROOT.TGraphErrors]] = []
  for ratioLabel, effInfosForLabel in effInfos.items():
    assert len(effInfosForLabel) == 2, f"Expect exactly 2 data samples to calculate ratio; but got {effInfosForLabel}"
    graphs: List[ROOT.TGraphErrors] = []
    for effInfo in effInfosForLabel.values():
      graphs.append(plotTools.getGraph1DFromValues(plotEfficiencies.getEffValuesForGraph1D(binningVar, effInfo)))
    ratioGraphs.append((ratioLabel, plotTools.calcRatioOfGraphs(graphs)))
  plotFitResults.plotGraphs1D(
    graphOrGraphs     = ratioGraphs,
    binningVar        = binningVar,
    yAxisTitle        = "Efficiency Ratio",
    pdfDirName        = pdfDirName,
    pdfFileBaseName   = "mm2_effratio",
    pdfFileNameSuffix = pdfFileNameSuffix,
    particle          = particle,
    channel           = channel,
    graphTitle        = f"{particle} Efficiency Ratio ({channel})" if graphTitle is None else graphTitle,
    graphMinimum      = 0.0,
    graphMaximum      = 1.3,
    skipBlack         = True if len(effInfos) > 1 else False,
    drawLegend        = True if len(effInfos) > 1 else False,
  )
  # draw line at 1
  #TODO implement line drawing routines in plotTools
  # xAxis = graph.GetXaxis()
  # yAxis = graph.GetYaxis()
  # oneLine = ROOT.TLine()
  # oneLine.SetLineStyle(ROOT.kDashed)
  # oneLine.SetLineColor(ROOT.kBlack)
  # oneLine.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 1, xAxis.GetBinUpEdge(xAxis.GetLast()), 1)


#TODO slice 2D graph instead
def overlayEfficiencyRatios2DSlices(
  effInfos:          Mapping[str, Mapping[Tuple[str, str], Sequence[EffInfo]]],
  binningVars:       Sequence[str],
  steppingVar:       str,
  pdfDirName:        str,
  pdfFileNameSuffix: str = "",
  particle:          str = "Proton",
  channel:           str = "4pi",
  graphTitle:        Optional[str] = None,
) -> None:
  """Plots efficiency ratios as a function of one binning variable while stepping through the bins of another variable given by `steppingVar` for all fits with matching 2D binning"""
  print(f"Plotting efficiency ratios for binning variables '{binningVars}' stepping through bins in '{steppingVar}'")
  ratioLabel, effInfosForLabel = next(iter(effInfos.items()))
  assert len(effInfosForLabel) == 2, f"Expect exactly 2 data samples to calculate ratio; but got {effInfosForLabel}"
  # filter efficiency infos that belong to given binning variables
  effInfos2D: Dict[Tuple[str, str], List[EffInfo]] = {}
  for key, efficiencies in effInfosForLabel.items():
    effInfos2D[key] = [
      efficiency for efficiency in efficiencies
      if (len(efficiency.binInfo.varNames) == 2) and (binningVars[0] in efficiency.binInfo.varNames) and (binningVars[1] in efficiency.binInfo.varNames)
    ]
    assert effInfos2D[key], f"Could not find any efficiencies that match binning variables '{binningVars}'"
  assert steppingVar in binningVars, f"Given stepping variable '{steppingVar}' must be in binning variables '{binningVars}'"
  assert steppingVar in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{steppingVar}'"
  xAxisVar = binningVars[0] if steppingVar == binningVars[1] else binningVars[1]
  firstKey = list(effInfos2D.keys())[0]
  # assume that binning info is identical for data sets
  steppingVarValues = set(effInfo.binInfo.centers[steppingVar] for effInfo in effInfos2D[firstKey])
  steppingVarWidth = effInfos2D[firstKey][0].binInfo.widths[steppingVar]  # assume equidistant binning
  for steppingVarValue in steppingVarValues:
    # construct input for 1D plotting function
    effInfos1D: Dict[Tuple[str, str], List[EffInfo]] = {}
    for key, efficiencies in effInfos2D.items():
      effInfos1D[key] = [
        EffInfo(BinInfo(
          name    = efficiency.binInfo.name,
          centers = {xAxisVar : efficiency.binInfo.centers[xAxisVar]},
          widths  = {xAxisVar : efficiency.binInfo.widths[xAxisVar]},
          dirName = efficiency.binInfo.dirName,
        ), efficiency.value) for efficiency in efficiencies if efficiency.binInfo.centers[steppingVar] == steppingVarValue
      ]
      assert effInfos1D[key], f"Could not find any efficiencies for bin center {steppingVar} == {steppingVarValue}"
    # plot current bin of stepping variable
    steppingRange = (f"{steppingVarValue - steppingVarWidth / 2.0}", f"{steppingVarValue + steppingVarWidth / 2.0}")
    steppingVarLabel = f"{steppingRange[0]} {BINNING_VAR_PLOT_INFO[steppingVar]['unit']} " \
      f"< {BINNING_VAR_PLOT_INFO[steppingVar]['label']} " \
      f"< {steppingRange[1]} {BINNING_VAR_PLOT_INFO[steppingVar]['unit']}"
    overlayEfficiencyRatios1D(
      effInfos          = {ratioLabel : effInfos1D},
      binningVar        = xAxisVar,
      pdfDirName        = pdfDirName,
      pdfFileNameSuffix = f"_{steppingVar}_{steppingRange[0]}_{steppingRange[1]}{pdfFileNameSuffix}",
      particle          = particle,
      channel           = channel,
      graphTitle        = f"{graphTitle}, {steppingVarLabel}",
    )


if __name__ == "__main__":
  makePlots.printGitInfo()
  ROOT.gROOT.SetBatch(True)
  plotTools.setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")

  ratiosToPlot: Dict[str, Tuple[str, str]] = {
    "2017_01-ver03" : (
      "./fits/2017_01-ver03/noShowers/BruFitOutput.bggen_2017_01-ver03_allFixed",
      "./fits/2017_01-ver03/noShowers/BruFitOutput.data_2017_01-ver03_allFixed",
    ),
    "2018_01-ver02" : (
      "./fits/2018_01-ver02/noShowers/BruFitOutput.bggen_2018_01-ver02_allFixed",
      "./fits/2018_01-ver02/noShowers/BruFitOutput.data_2018_01-ver02_allFixed",
    ),
    "2018_08-ver02" : (
      "./fits/2018_08-ver02/noShowers/BruFitOutput.bggen_2018_08-ver02_allFixed",
      "./fits/2018_08-ver02/noShowers/BruFitOutput.data_2018_08-ver02_allFixed",
    ),
    "2019_11-ver01" : (
      "./fits/2019_11-ver01/noShowers/BruFitOutput.bggen_2019_11-ver01_allFixed",
      "./fits/2019_11-ver01/noShowers/BruFitOutput.data_2019_11-ver01_allFixed",
    ),
  }

  pdfDirName = makePlots.makeDirPath("./ratios")
  graphTitle = "bggen MC / Real Data"
  effInfos:    Dict[str, Dict[Tuple[str, str], List[EffInfo]]] = {}
  binVarNames: Dict[str, Optional[List[Tuple[str, ...]]]]      = {}
  for ratioLabel, fitResults in ratiosToPlot.items():
    effInfos[ratioLabel], binVarNames[ratioLabel] = overlayEfficiencies.getEfficiencies(fitResultDirNames = tuple(fitResult for fitResult in fitResults))
  print("Plotting efficiency ratios")
  if effInfos and binVarNames:
    firstBinVarNames = next(iter(binVarNames.values()))  # get first entry
    for ratioLabel, binVarNamesForLabel in binVarNames.items():
      assert binVarNamesForLabel == firstBinVarNames, f"Data samples have different binnings: '{ratioLabel}' = {binVarNamesForLabel} vs. '{next(iter(binVarNames.keys()))}' = {firstBinVarNames}"
    if firstBinVarNames is not None:
      for binningVars in firstBinVarNames:
        if len(binningVars) == 1:
          overlayEfficiencyRatios1D(effInfos, binningVars[0], pdfDirName, graphTitle = graphTitle)
        if len(binningVars) == 2:
          overlayEfficiencyRatios2DSlices(effInfos, binningVars = binningVars[:2], steppingVar = binningVars[1], pdfDirName = pdfDirName, graphTitle = graphTitle)
