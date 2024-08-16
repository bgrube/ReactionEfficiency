#!/usr/bin/env python3


from __future__ import annotations

import argparse
from collections import defaultdict
from collections.abc import (
  Iterable,
  Mapping,
  Sequence,
)
from dataclasses import dataclass
import functools
import os
import sys

from uncertainties import UFloat, ufloat

import ROOT

from plotEfficiencies import (
  # BinInfo,
  BinningInfo,
  calculateEfficiencies,
  EffInfo,
  getEffValuesForGraph1D,
  getEffValuesForGraph2D,
  readYieldInfosForBinning,
)
from plotFitResults import (
  getAxisInfoForBinningVar,
  getBinningInfosFromDir,
  ParInfo,
  plotGraphs1D,
)
import plotTools
from plotTools import (
  getGraph1DFromValues,
  getGraph2DFromValues,
  makeDirPath,
  printGitInfo,
  setupPlotStyle,
  slice2DGraph,
)


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


def getEfficiencies(
  fitResultDirNames: Sequence[str],
  fitLabels:         Sequence[str] = [],
  dataSets:          Iterable[str] = ["Total", "Found", "Missing"],
  fitVariable:       str           = "MissingMassSquared_Measured",  # name of the fit variable
  useMissing:        bool          = True,  # True: 'Missing' dataset is used for efficiency calculation; False: 'Total' dataset is used instead
) -> tuple[dict[tuple[str, str], list[EffInfo]], list[BinningInfo | None]]:
  """Reads yields from given fit directories and calculates efficiencies for each directory"""
  assert len(fitResultDirNames) == len(set(fitResultDirNames)), f"List of fit-result directory names '{fitResultDirNames}' has duplicate elements"
  assert (not fitLabels) or len(fitLabels) == len(fitResultDirNames), f"Number of given fit labels ({len(fitLabels)}) does not match number of fit directories ({len(fitResultDirNames)})"
  print("Reading yields and calculating efficiencies")
  effInfos:     dict[tuple[str, str], list[EffInfo]] = {}    # effInfos[(<fit directory>, <fit label>)][<bin>]
  binningInfos: list[BinningInfo | None]             = []
  for fitIndex, fitResultDirName in enumerate(fitResultDirNames):
    yieldInfos: dict[str, list[ParInfo]] = {}  # yieldInfos[<dataset>][<bin>]
    for dataSet in dataSets:
      print(f"Reading yields for '{dataSet}' dataset")
      yieldInfos[dataSet]  = []
      binVarNamesInDataSet = []
      binningInfosForDataSet: list[BinningInfo | None] = getBinningInfosFromDir(f"{fitResultDirName}/{dataSet}")
      for binningInfo in binningInfosForDataSet:
        if binningInfo:
          yieldInfos[dataSet].extend(readYieldInfosForBinning(binningInfo, readIntegrals = False, fitVariable = fitVariable))  # append yield values
      if binningInfos:
        assert len(binningInfos) == len(binningInfosForDataSet), f"The number of binnings {len(binningInfosForDataSet)} for dataset '{dataSet}' is different from the number of binnings {len(binningInfos)} of the previous dataset"
        for binningInfo, binningInfoDataSet in zip(binningInfos, binningInfosForDataSet):
          assertMessage = f"The binning {binningInfoDataSet} for dataset '{dataSet}' is different from the binning {binningInfo} of the previous dataset"
          if binningInfo is None:
            assert binningInfoDataSet is None, assertMessage
          if binningInfoDataSet is None:
            assert binningInfo is None, assertMessage
          assert (binningInfo is not None and binningInfoDataSet is not None) and binningInfo.isSameBinningAs(binningInfoDataSet), assertMessage
      else:
        binningInfos = binningInfosForDataSet

    fitLabel = fitLabels[fitIndex] if fitLabels else fitResultDirName
    effInfos[(fitResultDirName, fitLabel)] = calculateEfficiencies(yieldInfos, useMissing)
  return effInfos, binningInfos


def overlayEfficiencies1D(
  effInfos:          Mapping[tuple[str, str], Sequence[EffInfo]],
  binningInfo:       BinningInfo,  # 1D binning information
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str        = "Proton_4pi_",
  pdfFileNameSuffix: str        = "",
  graphTitle:        str | None = None,
  skipBlack:         bool       = True,
) -> None:
  """Overlays efficiencies as a function of given binning variable for 1-dimensional binning"""
  assert len(binningInfo.varNames) == 1, f"Need 1-dimensional binning, but got {binningInfo}"
  binningVar = binningInfo.varNames[0]
  print(f"Overlaying efficiencies for binning variable '{binningVar}'")
  graphs1D: list[tuple[str, ROOT.TGraphErrors]] = [
    (fitLabel, getGraph1DFromValues(getEffValuesForGraph1D(binningVar, efficiencies)))
    for (_, fitLabel), efficiencies in effInfos.items()
  ]
  plotGraphs1D(
    graphOrGraphs     = graphs1D,
    binningInfo       = binningInfo,
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
  effInfos:          Mapping[tuple[str, str], Sequence[EffInfo]],
  binningInfo:       BinningInfo,  # 2D binning information
  steppingVar:       str,
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str  = "Proton_4pi_",
  pdfFileNameSuffix: str  = "",
  skipBlack:         bool = True,
) -> None:
  """Overlays efficiencies as a function of one binning variable while stepping through the bins of another variable given by `steppingVar` for all fits with matching 2D binning"""
  assert len(binningInfo.varNames) == 2, f"Need 2-dimensional binning, but got {binningInfo}"
  binningVars = binningInfo.varNames[:2]
  print(f"Overlaying efficiencies for binning variables '{binningVars}' stepping through bins in '{steppingVar}'")
  assert steppingVar in binningVars, f"Given stepping variable '{steppingVar}' must be in binning variables '{binningVars}'"
  _, steppingVarLabel, steppingVarUnit = getAxisInfoForBinningVar(steppingVar)
  binningVarIndex  = 0 if steppingVar == binningVars[1] else 1
  # read efficiency values and slice them into 1D graphs
  graphsToOverlay: defaultdict[tuple[float, float], list[tuple[str, ROOT.TGraphErrors]]] = defaultdict(list)
  for (_, fitLabel), efficiencies in effInfos.items():
    graphs1D: dict[tuple[float, float], ROOT.TGraphErrors] = slice2DGraph(
      getGraph2DFromValues(getEffValuesForGraph2D(binningVars, efficiencies)),
      plotTools.Graph2DVar.x if steppingVar == binningVars[0] else plotTools.Graph2DVar.y)
    for steppingVarBinRange, graph in graphs1D.items():
      graphsToOverlay[steppingVarBinRange].append((fitLabel, graph))
  for steppingVarBinRange, graphs in graphsToOverlay.items():
    # overlay 1D graphs for current bin of stepping variable
    steppingVarTitle = f"{steppingVarBinRange[0]} {steppingVarUnit} < {steppingVarLabel} " \
                       f"< {steppingVarBinRange[1]} {steppingVarUnit}"
    plotGraphs1D(
      graphOrGraphs     = graphs,
      binningInfo       = binningInfo.varBinningInfo(binningVars[binningVarIndex]),
      yAxisTitle        = "Efficiency",
      pdfDirName        = pdfDirName,
      pdfFileBaseName   = "mm2_eff",
      pdfFileNamePrefix = pdfFileNamePrefix,
      pdfFileNameSuffix = f"_{steppingVar}_{steppingVarBinRange[0]}_{steppingVarBinRange[1]}{pdfFileNameSuffix}",
      graphTitle        = steppingVarTitle,
      graphMinimum      = 0.0,
      graphMaximum      = 1.0,
      skipBlack         = skipBlack,
    )


if __name__ == "__main__":
  printGitInfo()
  ROOT.gROOT.SetBatch(True)
  setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Overlays efficiency graphs obtained from BruFit results in given directories using given labels.")
  parser.add_argument("--fitResult",         type = str, nargs = 2, metavar = ("DIR_PATH", "LABEL"), action = "append", help = "The path to the BruFit output directory of the fit that should be added to the overlay and the corresponding legend label; can be specified multiple times")
  parser.add_argument("--pdfFileNameSuffix", type = str, default = "", help = "PDF file-name suffix; (default: '%(default)s')")
  parser.add_argument("--useTotal",          action = "store_true", help = "If set, 'Total' distributions are used for efficiency calculation")
  args = parser.parse_args()

  fitRootDir = "./fits.nominal"
  pdfDirName = makeDirPath("./overlays")
  # fitRootDir = "./fits.pionComparison"
  # pdfDirName = makeDirPath("./overlays.pionComparison")
  skipBlack  = True
  # skipBlack  = False

  resultsToOverlay: dict[str, tuple[tuple[str, str], ...]] = {  # dict key is PDF file-name suffix
    # "bggen" : (
    #   (f"{fitRootDir}/2017_01-ver03/BruFitOutput.bggen_2017_01-ver03_allFixed", "2017_01-ver03"),
    #   (f"{fitRootDir}/2018_01-ver02/BruFitOutput.bggen_2018_01-ver02_allFixed", "2018_01-ver02"),
    #   (f"{fitRootDir}/2018_08-ver02/BruFitOutput.bggen_2018_08-ver02_allFixed", "2018_08-ver02"),
    #   (f"{fitRootDir}/2019_11-ver01/BruFitOutput.bggen_2019_11-ver01_allFixed", "2019_11-ver01"),
    # ),
    # "data" : (
    #   (f"{fitRootDir}/2017_01-ver03/BruFitOutput.data_2017_01-ver03_allFixed",  "2017_01-ver03"),
    #   (f"{fitRootDir}/2018_01-ver02/BruFitOutput.data_2018_01-ver02_allFixed",  "2018_01-ver02"),
    #   (f"{fitRootDir}/2018_08-ver02/BruFitOutput.data_2018_08-ver02_allFixed",  "2018_08-ver02"),
    #   (f"{fitRootDir}/2019_11-ver01/BruFitOutput.data_2019_11-ver01_allFixed",  "2019_11-ver01"),
    # ),
    "2017_01-ver03" : (
      (f"{fitRootDir}/2017_01-ver03/BruFitOutput.bggen_2017_01-ver03_allFixed", "bggen MC"),
      (f"{fitRootDir}/2017_01-ver03/BruFitOutput.data_2017_01-ver03_allFixed",  "Real Data"),
    ),
    # "2018_01-ver02" : (
    #   (f"{fitRootDir}/2018_01-ver02/BruFitOutput.bggen_2018_01-ver02_allFixed", "bggen MC"),
    #   (f"{fitRootDir}/2018_01-ver02/BruFitOutput.data_2018_01-ver02_allFixed",  "Real Data"),
    # ),
    # "2018_08-ver02" : (
    #   (f"{fitRootDir}/2018_08-ver02/BruFitOutput.bggen_2018_08-ver02_allFixed", "bggen MC"),
    #   (f"{fitRootDir}/2018_08-ver02/BruFitOutput.data_2018_08-ver02_allFixed",  "Real Data"),
    # ),
    # "2019_11-ver01" : (
    #   (f"{fitRootDir}/2019_11-ver01/BruFitOutput.bggen_2019_11-ver01_allFixed", "bggen MC"),
    #   (f"{fitRootDir}/2019_11-ver01/BruFitOutput.data_2019_11-ver01_allFixed",  "Real Data"),
    # ),
  }
  if args.fitResult:
    # read info from command-line arguments instead
    resultsToOverlay = {args.pdfFileNameSuffix : tuple((fitResult[0], fitResult[1]) for fitResult in args.fitResult)}

  for pdfFileNameSuffix, fitResults in resultsToOverlay.items():
    fitResultDirNames = tuple(fitResult[0] for fitResult in fitResults)
    fitLabels         = tuple(fitResult[1] for fitResult in fitResults)
    effInfos, binningInfos = getEfficiencies(fitResultDirNames, fitLabels, useMissing = not args.useTotal)
    print("Overlaying efficiencies")
    if effInfos and binningInfos:
      for binningInfo in binningInfos:
        if binningInfo:
          if len(binningInfo.varNames) == 1:
            overlayEfficiencies1D(effInfos, binningInfo, pdfDirName, f"_{pdfFileNameSuffix}", skipBlack = skipBlack)
          elif len(binningInfo.varNames) == 2:
            overlayEfficiencies2DSlices(effInfos, binningInfo, steppingVar = binningInfo.varNames[1],
              pdfDirName = pdfDirName, pdfFileNameSuffix = f"_{pdfFileNameSuffix}", skipBlack = skipBlack)
