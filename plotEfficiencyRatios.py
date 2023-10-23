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
from plotEfficiencies import EffInfo
from plotFitResults import BINNING_VAR_PLOT_INFO
import plotTools


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


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
  """Plots ratios of efficiencies as a function of `binningVar` for all given fits with 1D binning"""
  print(f"Plotting efficiency ratio for binning variable '{binningVar}'")
  assert len(effInfos) == 2, f"Expect exactly 2 data samples to calculate ratio: f{effInfos}"
  effValues: List[List[Tuple[UFloat, UFloat]]] = [[], []]
  for resultIndex, key in enumerate(effInfos.keys()):
    effValues[resultIndex] = (plotEfficiencies.getEffValuesForGraph1D(binningVar, effInfos[key]))
  # print(f"!!!0 {len(effValues[0])}: {effValues[0]}")
  # print(f"!!!1 {len(effValues[1])}: {effValues[1]}")
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
  # print(f"!!!ratios {len(ratios)}: {ratios}")
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
  graph.SetMaximum(1.5)
  # print(f"!!! {particle}_{channel}_mm2_effratio_{binningVar}_{ratioLabel}{pdfFileNameSuffix}")
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


if __name__ == "__main__":
  makePlots.printGitInfo()
  ROOT.gROOT.SetBatch(True)
  plotTools.setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")

  ratiosToPlot: Dict[str, Tuple[str, str]] = {
    # "2017_01-ver03" : (
    #   "./fits/2017_01-ver03/noShowers/BruFitOutput.bggen_2017_01-ver03_allFixed",
    #   "./fits/2017_01-ver03/noShowers/BruFitOutput.data_2017_01-ver03_allFixed",
    # ),
    # "2018_01-ver02" : (
    #   "./fits/2018_01-ver02/noShowers/BruFitOutput.bggen_2018_01-ver02_allFixed",
    #   "./fits/2018_01-ver02/noShowers/BruFitOutput.data_2018_01-ver02_allFixed",
    # ),
    # "2018_08-ver02" : (
    #   "./fits/2018_08-ver02/noShowers/BruFitOutput.bggen_2018_08-ver02_allFixed",
    #   "./fits/2018_08-ver02/noShowers/BruFitOutput.data_2018_08-ver02_allFixed",
    # ),
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
        # if len(binningVars) == 2:
        #   overlayEfficiencies2D(effInfos, binningVars = binningVars[:2], steppingVar = binningVars[1],
        #                         pdfDirName = pdfDirName, pdfFileNameSuffix = f"_{pdfFileNameSuffix}", skipBlack = skipBlack)
