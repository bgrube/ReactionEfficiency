#!/usr/bin/env python3


import argparse
from collections import defaultdict
from dataclasses import dataclass
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
  Set,
  Tuple,
)

from uncertainties import UFloat, ufloat

import ROOT

import makePlots
import plotEfficiencies
from plotEfficiencies import EffInfo
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


def overlayEfficiencies1D(
  effInfos:          Mapping[Tuple[str, str], Sequence[EffInfo]],
  binningVar:        str,
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str = "Proton_4pi_",
  pdfFileNameSuffix: str = "",
  graphTitle:        Optional[str] = None,
  skipBlack:         bool = True,
) -> None:
  """Overlays efficiencies as a function of `binningVar` for all given fits with 1D binning"""
  print(f"Overlaying efficiencies for binning variable '{binningVar}'")
  graphs1D: List[Tuple[str, ROOT.TGraphErrors]] = [(fitLabel, plotTools.getGraph1DFromValues(plotEfficiencies.getEffValuesForGraph1D(binningVar, efficiencies)))
                                                   for (_, fitLabel), efficiencies in effInfos.items()]
  plotFitResults.plotGraphs1D(
    graphs1D,
    binningVar,
    yAxisTitle        = "Efficiency",
    pdfDirName        = pdfDirName,
    pdfFileBaseName   = "mm2_eff",
    pdfFileNamePrefix = pdfFileNamePrefix,
    pdfFileNameSuffix = pdfFileNameSuffix,
    graphTitle        = f"Track-Finding Efficiency" if graphTitle is None else graphTitle,
    graphMinimum      = 0.0,
    graphMaximum      = 1.0,
    skipBlack         = skipBlack,
  )


def overlayEfficiencies2DSlices(
  effInfos:          Mapping[Tuple[str, str], Sequence[EffInfo]],
  binningVars:       Sequence[str],
  steppingVar:       str,
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str = "Proton_4pi_",
  pdfFileNameSuffix: str = "",
  skipBlack:         bool = True,
) -> None:
  """Overlays efficiencies as a function of one binning variable while stepping through the bins of another variable given by `steppingVar` for all fits with matching 2D binning"""
  print(f"Overlaying efficiencies for binning variables '{binningVars}' stepping through bins in '{steppingVar}'")
  assert steppingVar in binningVars, f"Given stepping variable '{steppingVar}' must be in binning variables '{binningVars}'"
  assert steppingVar in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{steppingVar}'"
  binningVarIndex  = 0 if steppingVar == binningVars[1] else 1
  # read efficiency values and slice them into 1D graphs
  graphsToOverlay: Dict[Tuple[float, float], List[Tuple[str, ROOT.TGraphErrors]]] = defaultdict(list)
  for (_, fitLabel), efficiencies in effInfos.items():
    graphs1D: Dict[Tuple[float, float], ROOT.TGraphErrors] = plotTools.slice2DGraph(
      plotTools.getGraph2DFromValues(plotEfficiencies.getEffValuesForGraph2D(binningVars, efficiencies)),
      plotTools.Graph2DVar.x if steppingVar == binningVars[0] else plotTools.Graph2DVar.y)
    for steppingVarBinRange, graph in graphs1D.items():
      graphsToOverlay[steppingVarBinRange].append((fitLabel, graph))
  for steppingVarBinRange, graphs in graphsToOverlay.items():
    # overlay 1D graphs for current bin of stepping variable
    steppingVarLabel = f"{steppingVarBinRange[0]} {BINNING_VAR_PLOT_INFO[steppingVar]['unit']} " \
      f"< {BINNING_VAR_PLOT_INFO[steppingVar]['label']} " \
      f"< {steppingVarBinRange[1]} {BINNING_VAR_PLOT_INFO[steppingVar]['unit']}"
    plotFitResults.plotGraphs1D(
      graphs,
      binningVars[binningVarIndex],
      yAxisTitle        = "Efficiency",
      pdfDirName        = pdfDirName,
      pdfFileBaseName   = "mm2_eff",
      pdfFileNamePrefix = pdfFileNamePrefix,
      pdfFileNameSuffix = f"_{steppingVar}_{steppingVarBinRange[0]}_{steppingVarBinRange[1]}{pdfFileNameSuffix}",
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
    "bggen" : (
      ("./fits/2017_01-ver03/noShowers/BruFitOutput.bggen_2017_01-ver03_allFixed", "2017_01-ver03"),
      ("./fits/2018_01-ver02/noShowers/BruFitOutput.bggen_2018_01-ver02_allFixed", "2018_01-ver02"),
      ("./fits/2018_08-ver02/noShowers/BruFitOutput.bggen_2018_08-ver02_allFixed", "2018_08-ver02"),
      ("./fits/2019_11-ver01/noShowers/BruFitOutput.bggen_2019_11-ver01_allFixed", "2019_11-ver01"),
    ),
    "data" : (
      ("./fits/2017_01-ver03/noShowers/BruFitOutput.data_2017_01-ver03_allFixed",  "2017_01-ver03"),
      ("./fits/2018_01-ver02/noShowers/BruFitOutput.data_2018_01-ver02_allFixed",  "2018_01-ver02"),
      ("./fits/2018_08-ver02/noShowers/BruFitOutput.data_2018_08-ver02_allFixed",  "2018_08-ver02"),
      ("./fits/2019_11-ver01/noShowers/BruFitOutput.data_2019_11-ver01_allFixed",  "2019_11-ver01"),
    ),
    "2017_01-ver03" : (
      ("./fits/2017_01-ver03/noShowers/BruFitOutput.bggen_2017_01-ver03_allFixed", "bggen MC"),
      ("./fits/2017_01-ver03/noShowers/BruFitOutput.data_2017_01-ver03_allFixed",  "Real Data"),
    ),
    "2018_01-ver02" : (
      ("./fits/2018_01-ver02/noShowers/BruFitOutput.bggen_2018_01-ver02_allFixed", "bggen MC"),
      ("./fits/2018_01-ver02/noShowers/BruFitOutput.data_2018_01-ver02_allFixed",  "Real Data"),
    ),
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
          overlayEfficiencies1D(effInfos, binningVars[0], pdfDirName, f"_{pdfFileNameSuffix}", skipBlack = skipBlack)
        if len(binningVars) == 2:
          overlayEfficiencies2DSlices(effInfos, binningVars = binningVars[:2], steppingVar = binningVars[1],
                                pdfDirName = pdfDirName, pdfFileNameSuffix = f"_{pdfFileNameSuffix}", skipBlack = skipBlack)
