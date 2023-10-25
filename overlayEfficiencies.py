#!/usr/bin/env python3


import argparse
from dataclasses import dataclass  # builtin in Python 3.7+
import functools
import numpy as np
import os
import sys
from typing import (
  Dict,
  Iterable,
  List,
  Mapping,
  Optional,
  Sequence,
  Tuple,
)

from uncertainties import UFloat, ufloat

import ROOT

import makePlots
import plotEfficiencies
from plotEfficiencies import EffInfo, BinInfo
import plotFitResults
from plotFitResults import ParInfo, BINNING_VAR_PLOT_INFO
import plotTools


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


def getEfficiencies(
  fitResultDirNames: Sequence[str],
  fitLabels:         Sequence[str] = [],
  dataSets:          Iterable[str] = ["Total", "Found", "Missing"],
) -> Tuple[Dict[Tuple[str, str], List[EffInfo]], Optional[List[Tuple[str, ...]]]]:
  """Reads yields from given fit directories and calculates efficiencies for each directory"""
  assert len(fitResultDirNames) == len(set(fitResultDirNames)), f"List of fit-result directory names '{fitResultDirNames}' has duplicate elements"
  assert (not fitLabels) or len(fitLabels) == len(fitResultDirNames), f"Number of given fit labels ({len(fitLabels)}) does not match number of fit directories ({len(fitResultDirNames)})"
  print("Reading yields and calculating efficiencies")
  effInfos:    Dict[Tuple[str, str], List[EffInfo]] = {}    # effInfos[(<fit directory>, <fit label>)][<bin>]
  binVarNames: Optional[List[Tuple[str, ...]]]      = None  # binning variables for each binning
  for fitIndex, fitResultDirName in enumerate(fitResultDirNames):
    yieldInfos: Dict[str, List[ParInfo]] = {}  # yieldInfos[<dataset>][<bin>]
    for dataSet in dataSets:
      print(f"Reading yields for '{dataSet}' dataset")
      yieldInfos[dataSet]  = []
      binVarNamesInDataSet = []
      for binningInfo in plotFitResults.getBinningInfosFromDir(f"{fitResultDirName}/{dataSet}"):
        if binningInfo:
          binVarNamesInDataSet.append(binningInfo.varNames)
          yieldInfos[dataSet][len(yieldInfos[dataSet]):] = plotEfficiencies.readYieldInfosForBinning(binningInfo)  # append yield values
      if binVarNames is None:
        binVarNames = binVarNamesInDataSet
      else:
        assert binVarNamesInDataSet == binVarNames, f"The binning variables {binVarNamesInDataSet} of dataset '{dataSet}' are different from the binning variables {binVarNames} of the previous one"
    fitLabel = fitLabels[fitIndex] if fitLabels else fitResultDirName
    effInfos[(fitResultDirName, fitLabel)] = plotEfficiencies.calculateEfficiencies(yieldInfos)
  return effInfos, binVarNames


def plotGraphs1D(
  graphs:            Sequence[Tuple[str, ROOT.TGraph]],
  binningVar:        str,
  yAxisTitle:        str,
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNameSuffix: str = "",
  particle:          str = "Proton",
  channel:           str = "4pi",
  graphTitle:        Optional[str] = None,
  graphMinimum:      float = 0.0,
  graphMaximum:      Optional[float] = None,
  skipBlack:         bool = True,
):
  """Overlays all given graphs"""
  multiGraph = ROOT.TMultiGraph()
  for styleIndex, (legendLabel, graph) in enumerate(graphs):
    graph.SetTitle(legendLabel)
    plotTools.setCbFriendlyStyle(graph, styleIndex, skipBlack = False if len(graphs) == 1 else skipBlack)
    multiGraph.Add(graph)
  if graphTitle is not None:
    multiGraph.SetTitle(graphTitle)  # !Note! if this is executed after setting axis titles, no title is printed; seems like a ROOT bug
  assert binningVar in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVar}'"
  multiGraph.GetXaxis().SetTitle(f"{BINNING_VAR_PLOT_INFO[binningVar]['label']} ({BINNING_VAR_PLOT_INFO[binningVar]['unit']})")
  multiGraph.GetYaxis().SetTitle(yAxisTitle)
  if graphMinimum is not None:
    multiGraph.SetMinimum(graphMinimum)
  if graphMaximum is not None:
    multiGraph.SetMaximum(graphMaximum)
  legendLabels = tuple(legendLabel.replace(' ', '_') for legendLabel, _ in graphs)  # reformat legend labels so that they can be used in PDF file name
  canv = ROOT.TCanvas(f"{particle}_{channel}_mm2_eff_{binningVar}_{'_'.join(legendLabels)}{pdfFileNameSuffix}", "")
  multiGraph.Draw("APZ")
  canv.BuildLegend()
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


def overlayEfficiencies1D(
  effInfos:          Mapping[Tuple[str, str], Sequence[EffInfo]],
  binningVar:        str,
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNameSuffix: str = "",
  particle:          str = "Proton",
  channel:           str = "4pi",
  graphTitle:        Optional[str] = None,
  skipBlack:         bool = True,
):
  #TODO use overlayGraphs1D
  """Overlays efficiencies as a function of `binningVar` for all given fits with 1D binning"""
  print(f"Overlaying efficiencies for binning variable '{binningVar}'")
  efficiencyMultiGraph = ROOT.TMultiGraph()
  efficiencyGraphs = {}  # store graphs here to keep them in memory
  shiftFraction = 0.0
  styleIndex = 0
  for (fitResultDirName, fitLabel), efficiencies in effInfos.items():
    graph = efficiencyGraphs[fitResultDirName] = plotFitResults.getParValueGraph1D(plotEfficiencies.getEffValuesForGraph1D(binningVar, efficiencies), shiftFraction)
    # shiftFraction += 0.01
    graph.SetTitle(fitLabel)
    plotTools.setCbFriendlyStyle(graph, styleIndex, skipBlack = False if len(effInfos) == 1 else skipBlack)
    styleIndex += 1
    efficiencyMultiGraph.Add(graph)
  if graphTitle is None:
    efficiencyMultiGraph.SetTitle(f"{particle} Track-Finding Efficiency ({channel})")
  else:
    efficiencyMultiGraph.SetTitle(graphTitle)
  assert binningVar in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVar}'"
  efficiencyMultiGraph.GetXaxis().SetTitle(f"{BINNING_VAR_PLOT_INFO[binningVar]['label']} ({BINNING_VAR_PLOT_INFO[binningVar]['unit']})")
  efficiencyMultiGraph.GetYaxis().SetTitle("Efficiency")
  efficiencyMultiGraph.SetMinimum(0)
  efficiencyMultiGraph.SetMaximum(1)
  fitLabels = tuple(key[1].replace(' ', '_') for key in effInfos.keys())
  canv = ROOT.TCanvas(f"{particle}_{channel}_mm2_eff_{binningVar}_{'_'.join(fitLabels)}{pdfFileNameSuffix}", "")
  efficiencyMultiGraph.Draw("APZ")
  canv.BuildLegend()
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


def overlayEfficiencies1DNew(
  effInfos:          Mapping[Tuple[str, str], Sequence[EffInfo]],
  binningVar:        str,
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNameSuffix: str = "",
  particle:          str = "Proton",
  channel:           str = "4pi",
  graphTitle:        Optional[str] = None,
  skipBlack:         bool = True,
):
  """Overlays efficiencies as a function of `binningVar` for all given fits with 1D binning"""
  print(f"Overlaying efficiencies for binning variable '{binningVar}'")
  graphs1D: List[Tuple[str, ROOT.TGraphErrors]] = []
  for (_, fitLabel), efficiencies in effInfos.items():
    graphs1D.append((fitLabel, plotFitResults.getParValueGraph1D(plotEfficiencies.getEffValuesForGraph1D(binningVar, efficiencies))))
  plotGraphs1D(
    graphs1D,
    binningVar,
    yAxisTitle        = "Efficiency",
    pdfDirName        = pdfDirName,
    pdfFileNameSuffix = pdfFileNameSuffix,
    particle          = particle,
    channel           = channel,
    graphTitle        = f"{particle} Track-Finding Efficiency ({channel})" if graphTitle is None else graphTitle,
    graphMinimum      = 0.0,
    graphMaximum      = 1.0,
    skipBlack         = skipBlack,
  )


def overlayEfficiencies2D(
  effInfos:          Mapping[Tuple[str, str], Sequence[EffInfo]],
  binningVars:       Sequence[str],
  steppingVar:       str,
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNameSuffix: str = "",
  particle:          str = "Proton",
  channel:           str = "4pi",
  skipBlack:         bool = True,
):
  """Overlays efficiencies as a function of one binning variable while stepping through the bins of another variable given by `steppingVar` for all fits with matching 2D binning"""
  print(f"Overlaying efficiencies for binning variables '{binningVars}' stepping through bins in '{steppingVar}'")
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
    overlayEfficiencies1D(
      effInfos = effInfos1D,
      binningVar = xAxisVar,
      pdfDirName = pdfDirName,
      pdfFileNameSuffix = f"_{steppingVar}_{steppingRange[0]}_{steppingRange[1]}{pdfFileNameSuffix}",
      particle = particle,
      channel = channel,
      graphTitle = steppingVarLabel,
      skipBlack = skipBlack,
    )


def overlayEfficiencies2DNew(
  effInfos:          Mapping[Tuple[str, str], Sequence[EffInfo]],
  binningVars:       Sequence[str],
  steppingVar:       str,
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNameSuffix: str = "",
  particle:          str = "Proton",
  channel:           str = "4pi",
  skipBlack:         bool = True,
):
  """Overlays efficiencies as a function of one binning variable while stepping through the bins of another variable given by `steppingVar` for all fits with matching 2D binning"""
  print(f"Overlaying efficiencies for binning variables '{binningVars}' stepping through bins in '{steppingVar}'")
  assert steppingVar in binningVars, f"Given stepping variable '{steppingVar}' must be in binning variables '{binningVars}'"
  assert steppingVar in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{steppingVar}'"
  steppingVarIndex = 0 if steppingVar == binningVars[0] else 1
  binningVarIndex  = 0 if steppingVarIndex == 1 else 1
  # read efficiency values assuming equidistant binning
  steppingVarBinCenters = set()
  steppingVarBinWidths  = set()
  graphValues: Dict[Tuple[str, str], Tuple[Tuple[UFloat, UFloat, UFloat], ...]] = {}  # [(dir name , label)][value index][0 = x, 1 = y, 2 = eff]
  for key, efficiencies in effInfos.items():
    values: Tuple[Tuple[UFloat, UFloat, UFloat], ...] = plotEfficiencies.getEffValuesForGraph2D(binningVars, efficiencies)
    if not values:  # nothing to plot
      continue
    steppingVarBinCenters.update((value[steppingVarIndex].nominal_value for value in values))
    steppingVarBinWidths.update ((value[steppingVarIndex].std_dev * 2   for value in values))
    graphValues[key] = values
  assert len(steppingVarBinWidths) == 1; f"Binning for stepping variable is not equidistant; found bin widths {steppingVarBinWidths}"
  steppingVarBinWidth = list(steppingVarBinWidths)[0]
  for steppingVarBinCenter in sorted(steppingVarBinCenters):
    # construct 1D graphs for current bin of stepping variable
    graphs1D: List[Tuple[str, ROOT.TGraphErrors]] = []
    for (fitResultDirName, fitLabel), values in graphValues.items():
      xVals = np.array([value[binningVarIndex].nominal_value for value in values if value[steppingVarIndex].nominal_value == steppingVarBinCenter], dtype = "d")
      xErrs = np.array([value[binningVarIndex].std_dev       for value in values if value[steppingVarIndex].nominal_value == steppingVarBinCenter], dtype = "d")
      yVals = np.array([value[2].nominal_value               for value in values if value[steppingVarIndex].nominal_value == steppingVarBinCenter], dtype = "d")
      yErrs = np.array([value[2].std_dev                     for value in values if value[steppingVarIndex].nominal_value == steppingVarBinCenter], dtype = "d")
      assert len(xVals) > 0, f"Could not find any efficiencies for '{(fitResultDirName, fitLabel)}' for bin center {steppingVar} == {steppingVarBinCenter}"
      graphs1D.append((fitLabel, ROOT.TGraphErrors(len(xVals), xVals, yVals, xErrs, yErrs)))
    # overlay 1D graphs for current bin of stepping variable
    steppingRange = (f"{steppingVarBinCenter - steppingVarBinWidth / 2.0}", f"{steppingVarBinCenter + steppingVarBinWidth / 2.0}")
    steppingVarLabel = f"{steppingRange[0]} {BINNING_VAR_PLOT_INFO[steppingVar]['unit']} " \
      f"< {BINNING_VAR_PLOT_INFO[steppingVar]['label']} " \
      f"< {steppingRange[1]} {BINNING_VAR_PLOT_INFO[steppingVar]['unit']}"
    plotGraphs1D(
      graphs1D,
      binningVars[binningVarIndex],
      yAxisTitle        = "Efficiency",
      pdfDirName        = pdfDirName,
      pdfFileNameSuffix = f"_{steppingVar}_{steppingRange[0]}_{steppingRange[1]}{pdfFileNameSuffix}",
      particle          = particle,
      channel           = channel,
      graphTitle        = steppingVarLabel,
      graphMinimum      = 0.0,
      graphMaximum      = 1.0,
      skipBlack         = skipBlack,
    )


if __name__ == "__main__":
  makePlots.printGitInfo()
  ROOT.gROOT.SetBatch(True)
  plotTools.setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Overlays efficiency graphs obtained from BruFit results in given directories using given labels.")
  parser.add_argument("--fitResult",         type = str, nargs = 2, metavar = ("DIR_PATH", "LABEL"), action = "append", help = "The path to the BruFit output directory of the fit that should be added to the overlay and the corresponding legend label; can be specified multiple times")
  parser.add_argument("--pdfFileNameSuffix", type = str, default = "", help = "PDF file-name suffix; (default: '%(default)s')")
  args = parser.parse_args()

  resultsToOverlay: Dict[str, Tuple[Tuple[str, str], ...]] = {  # dict key is PDF file-name suffix
    # "bggen" : (
    #   ("./fits/2017_01-ver03/noShowers/BruFitOutput.bggen_2017_01-ver03_allFixed", "2017_01-ver03"),
    #   ("./fits/2018_01-ver02/noShowers/BruFitOutput.bggen_2018_01-ver02_allFixed", "2018_01-ver02"),
    #   ("./fits/2018_08-ver02/noShowers/BruFitOutput.bggen_2018_08-ver02_allFixed", "2018_08-ver02"),
    #   ("./fits/2019_11-ver01/noShowers/BruFitOutput.bggen_2019_11-ver01_allFixed", "2019_11-ver01"),
    # ),
    # "data" : (
    #   ("./fits/2017_01-ver03/noShowers/BruFitOutput.data_2017_01-ver03_allFixed",  "2017_01-ver03"),
    #   ("./fits/2018_01-ver02/noShowers/BruFitOutput.data_2018_01-ver02_allFixed",  "2018_01-ver02"),
    #   ("./fits/2018_08-ver02/noShowers/BruFitOutput.data_2018_08-ver02_allFixed",  "2018_08-ver02"),
    #   ("./fits/2019_11-ver01/noShowers/BruFitOutput.data_2019_11-ver01_allFixed",  "2019_11-ver01"),
    # ),
    # "2017_01-ver03" : (
    #   ("./fits/2017_01-ver03/noShowers/BruFitOutput.bggen_2017_01-ver03_allFixed", "bggen MC"),
    #   ("./fits/2017_01-ver03/noShowers/BruFitOutput.data_2017_01-ver03_allFixed",  "Real Data"),
    # ),
    # "2018_01-ver02" : (
    #   ("./fits/2018_01-ver02/noShowers/BruFitOutput.bggen_2018_01-ver02_allFixed", "bggen MC"),
    #   ("./fits/2018_01-ver02/noShowers/BruFitOutput.data_2018_01-ver02_allFixed",  "Real Data"),
    # ),
    "2018_08-ver02" : (
      ("./fits/2018_08-ver02/noShowers/BruFitOutput.bggen_2018_08-ver02_allFixed", "bggen MC"),
      ("./fits/2018_08-ver02/noShowers/BruFitOutput.data_2018_08-ver02_allFixed",  "Real Data"),
    ),
    "2019_11-ver01" : (
      ("./fits/2019_11-ver01/noShowers/BruFitOutput.bggen_2019_11-ver01_allFixed", "bggen MC"),
      ("./fits/2019_11-ver01/noShowers/BruFitOutput.data_2019_11-ver01_allFixed",  "Real Data"),
    ),
  }
  if args.fitResult:
    # read info from command-line arguments instead
    resultsToOverlay = {args.pdfFileNameSuffix : tuple((fitResult[0], fitResult[1]) for fitResult in args.fitResult)}

  skipBlack = True
  # skipBlack = False
  pdfDirName = makePlots.makeDirPath("./overlays")
  for pdfFileNameSuffix, fitResults in resultsToOverlay.items():
    fitResultDirNames = tuple(fitResult[0] for fitResult in fitResults)
    fitLabels         = tuple(fitResult[1] for fitResult in fitResults)
    effInfos, binVarNames = getEfficiencies(fitResultDirNames, fitLabels)
    print("Overlaying efficiencies")
    if effInfos and binVarNames:
      for binningVars in binVarNames:
        if len(binningVars) == 1:
          # overlayEfficiencies1D(effInfos, binningVars[0], pdfDirName, f"_{pdfFileNameSuffix}", skipBlack = skipBlack)
          overlayEfficiencies1DNew(effInfos, binningVars[0], pdfDirName, f"_{pdfFileNameSuffix}", skipBlack = skipBlack)
        if len(binningVars) == 2:
          # overlayEfficiencies2D(effInfos, binningVars = binningVars[:2], steppingVar = binningVars[1],
          #                       pdfDirName = pdfDirName, pdfFileNameSuffix = f"_{pdfFileNameSuffix}", skipBlack = skipBlack)
          overlayEfficiencies2DNew(effInfos, binningVars = binningVars[:2], steppingVar = binningVars[1],
                                pdfDirName = pdfDirName, pdfFileNameSuffix = f"_{pdfFileNameSuffix}", skipBlack = skipBlack)
