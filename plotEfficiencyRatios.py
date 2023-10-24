#!/usr/bin/env python3


import functools
import numpy as np
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
from plotFitResults import BINNING_VAR_PLOT_INFO
import plotTools


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


#TODO simplify by just creating ratio graph and using existing plot function
def plotEfficiencyRatio1D(
  effInfos:          Mapping[Tuple[str, str], Sequence[EffInfo]],
  binningVar:        str,
  ratioLabel:        str,
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNameSuffix: str = "",
  particle:          str = "Proton",
  channel:           str = "4pi",
  graphTitle:        Optional[str] = None,
):
  """Plots efficiency ratios as a function of `binningVar` for all given fits with 1D binning"""
  print(f"Plotting efficiency ratio for binning variable '{binningVar}'")
  assert len(effInfos) == 2, f"Expect exactly 2 data samples to calculate ratio: f{effInfos}"
  effValues: List[List[Tuple[UFloat, UFloat]]] = [[], []]
  for resultIndex, key in enumerate(effInfos.keys()):
    effValues[resultIndex] = (plotEfficiencies.getEffValuesForGraph1D(binningVar, effInfos[key]))
  ratios: List[Tuple[UFloat, UFloat]] = []
  for effIndex in range(len(effValues[0])):
    binVal  = effValues[0][effIndex][0]
    effVal1 = effValues[0][effIndex][1]
    # it is not guaranteed that the two data sets contain efficiencies for all bins
    # find corresponding bin value in other data set
    effVal2 = None
    for effValue in effValues[1]:
      if effValue[0].nominal_value == binVal.nominal_value and effValue[0].std_dev == binVal.std_dev:
        effVal2 = effValue[1]
        break
    if effVal2 is None:
      continue
    # calculate ratio
    ratios.append((binVal, effVal1 / effVal2))
  # create graph
  xVals = np.array([ratio[0].nominal_value for ratio in ratios], dtype = "d")
  xErrs = np.array([ratio[0].std_dev       for ratio in ratios], dtype = "d")
  yVals = np.array([ratio[1].nominal_value for ratio in ratios], dtype = "d")
  yErrs = np.array([ratio[1].std_dev       for ratio in ratios], dtype = "d")
  graph = ROOT.TGraphErrors(len(xVals), xVals, yVals, xErrs, yErrs)
  plotTools.setCbFriendlyStyle(graph, 0, skipBlack = False)
  if graphTitle is None:
    graph.SetTitle(f"{particle} Efficiency Ratio ({channel})")
  else:
    graph.SetTitle(graphTitle)
  assert binningVar in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVar}'"
  graph.GetXaxis().SetTitle(f"{BINNING_VAR_PLOT_INFO[binningVar]['label']} ({BINNING_VAR_PLOT_INFO[binningVar]['unit']})")
  graph.GetYaxis().SetTitle("Efficiency Ratio")
  graph.SetMinimum(0)
  graph.SetMaximum(1.3)
  canv = ROOT.TCanvas(f"{particle}_{channel}_mm2_effratio_{binningVar}_{ratioLabel}{pdfFileNameSuffix}", "")
  graph.Draw("APZ")
  # draw line at 1
  #TODO implement line drawing routines in plotTools
  xAxis = graph.GetXaxis()
  yAxis = graph.GetYaxis()
  oneLine = ROOT.TLine()
  oneLine.SetLineStyle(ROOT.kDashed)
  oneLine.SetLineColor(ROOT.kBlack)
  oneLine.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 1, xAxis.GetBinUpEdge(xAxis.GetLast()), 1)
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


#TODO slice 2D graph instead
def plotEfficiencyRatio2D(
  effInfos:          Mapping[Tuple[str, str], Sequence[EffInfo]],
  binningVars:       Sequence[str],
  steppingVar:       str,
  ratioLabel:        str,
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNameSuffix: str = "",
  particle:          str = "Proton",
  channel:           str = "4pi",
  graphTitle:        Optional[str] = None,
):
  """Plots efficiency ratios as a function of one binning variable while stepping through the bins of another variable given by `steppingVar` for all fits with matching 2D binning"""
  print(f"Plotting efficiency ratios for binning variables '{binningVars}' stepping through bins in '{steppingVar}'")
  # filter efficiency infos that belong to given binning variables
  effInfos2D: Dict[Tuple[str, str], List[EffInfo]] = {}
  for key, efficiencies in effInfos.items():
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
    plotEfficiencyRatio1D(
      effInfos          = effInfos1D,
      binningVar        = xAxisVar,
      ratioLabel        = ratioLabel,
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
  for ratioLabel, fitResults in ratiosToPlot.items():
    effInfos, binVarNames = overlayEfficiencies.getEfficiencies(fitResultDirNames = tuple(fitResult for fitResult in fitResults))
    print("Plotting efficiency ratios")
    if effInfos and binVarNames:
      for binningVars in binVarNames:
        if len(binningVars) == 1:
          plotEfficiencyRatio1D(effInfos, binningVars[0], ratioLabel, pdfDirName, graphTitle = graphTitle)
        if len(binningVars) == 2:
          plotEfficiencyRatio2D(effInfos, binningVars = binningVars[:2], steppingVar = binningVars[1],
                                ratioLabel = ratioLabel, pdfDirName = pdfDirName, graphTitle = graphTitle)
