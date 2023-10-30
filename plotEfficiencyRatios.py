#!/usr/bin/env python3


from collections import defaultdict
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

import overlayEfficiencies
import plotEfficiencies
from plotEfficiencies import BinInfo, EffInfo
import plotFitResults
from plotFitResults import BINNING_VAR_PLOT_INFO
import plotTools


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


def overlayEfficiencyRatios1D(
  effInfos:          Mapping[str, Mapping[Tuple[str, str], Sequence[EffInfo]]],  # [ratioLabel][(fitResultDirName, fitLabel)][bin index]
  binningVar:        str,
  pdfDirName:        str,
  pdfFileNamePrefix: str = "Proton_4pi_",
  pdfFileNameSuffix: str = "",
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
    ratioGraphs.append((ratioLabel, plotTools.calcRatioOfGraphs1D(graphs)))
  plotFitResults.plotGraphs1D(
    graphOrGraphs     = ratioGraphs,
    binningVar        = binningVar,
    yAxisTitle        = "Efficiency Ratio",
    pdfDirName        = pdfDirName,
    pdfFileBaseName   = "mm2_effratio",
    pdfFileNamePrefix = pdfFileNamePrefix,
    pdfFileNameSuffix = pdfFileNameSuffix,
    graphTitle        = f"Efficiency Ratio" if graphTitle is None else graphTitle,
    graphMinimum      = 0.0,
    graphMaximum      = 1.3,
    skipBlack         = True if len(ratioGraphs) > 1 else False,
    drawLegend        = True if len(ratioGraphs) > 1 else False,
  )
  # draw line at 1
  #TODO implement line drawing routines in plotTools
  # xAxis = graph.GetXaxis()
  # yAxis = graph.GetYaxis()
  # oneLine = ROOT.TLine()
  # oneLine.SetLineStyle(ROOT.kDashed)
  # oneLine.SetLineColor(ROOT.kBlack)
  # oneLine.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 1, xAxis.GetBinUpEdge(xAxis.GetLast()), 1)

def overlayEfficiencyRatios2DSlices(
  effInfos:          Mapping[str, Mapping[Tuple[str, str], Sequence[EffInfo]]],  # [ratioLabel][(fitResultDirName, fitLabel)][bin index]
  binningVars:       Sequence[str],
  steppingVar:       str,
  pdfDirName:        str,
  pdfFileNamePrefix: str = "Proton_4pi_",
  pdfFileNameSuffix: str = "",
  graphTitle:        Optional[str] = None,
) -> None:
  """Plots efficiency ratios as a function of one binning variable while stepping through the bins of another variable given by `steppingVar` for all fits with matching 2D binning"""
  print(f"Plotting efficiency ratios for binning variables '{binningVars}' stepping through bins in '{steppingVar}'")
  assert steppingVar in binningVars, f"Given stepping variable '{steppingVar}' must be in binning variables '{binningVars}'"
  assert steppingVar in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{steppingVar}'"
  binningVarIndex  = 0 if steppingVar == binningVars[1] else 1
  graphsToOverlay: Dict[Tuple[float, float], List[Tuple[str, ROOT.TGraphErrors]]] = defaultdict(list)
  for ratioLabel, effInfosForLabel in effInfos.items():
    assert len(effInfosForLabel) == 2, f"Expect exactly 2 data samples to calculate ratio; but got {effInfosForLabel}"
    # get efficiencies as 2D graphs
    efficiencyGraphs2D: Tuple[ROOT.TGraph2DErrors, ...] = tuple(plotTools.getGraph2DFromValues(plotEfficiencies.getEffValuesForGraph2D(binningVars, efficiencies))
                                                                for efficiencies in effInfosForLabel.values())
    # calculate 2D graph for efficiency ratios and slice it to 1D graphs
    ratioGraphs1D: Dict[Tuple[float, float], ROOT.TGraphErrors] = plotTools.slice2DGraph(plotTools.calcRatioOfGraphs2D(efficiencyGraphs2D),
      plotTools.Graph2DVar.x if steppingVar == binningVars[0] else plotTools.Graph2DVar.y)
    for steppingVarBinRange, graph in ratioGraphs1D.items():
      graphsToOverlay[steppingVarBinRange].append((ratioLabel, graph))
  for steppingVarBinRange, graphs in graphsToOverlay.items():
    # overlay 1D graphs for current bin of stepping variable
    steppingVarLabel = f"{steppingVarBinRange[0]} {BINNING_VAR_PLOT_INFO[steppingVar]['unit']} " \
      f"< {BINNING_VAR_PLOT_INFO[steppingVar]['label']} " \
      f"< {steppingVarBinRange[1]} {BINNING_VAR_PLOT_INFO[steppingVar]['unit']}"
    plotFitResults.plotGraphs1D(
      graphs,
      binningVars[binningVarIndex],
      yAxisTitle        = "Efficiency Ratio",
      pdfDirName        = pdfDirName,
      pdfFileBaseName   = "mm2_effratio",
      pdfFileNamePrefix = pdfFileNamePrefix,
      pdfFileNameSuffix = f"_{steppingVar}_{steppingVarBinRange[0]}_{steppingVarBinRange[1]}{pdfFileNameSuffix}",
      graphTitle        = f"{graphTitle}, {steppingVarLabel}",
      graphMinimum      = 0.0,
      graphMaximum      = 1.3,
      skipBlack         = True if len(graphs) > 1 else False,
      drawLegend        = True if len(graphs) > 1 else False,
    )


if __name__ == "__main__":
  plotTools.printGitInfo()
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

  pdfDirName = plotTools.makeDirPath("./ratios")
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
