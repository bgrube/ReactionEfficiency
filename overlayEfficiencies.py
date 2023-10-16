#!/usr/bin/env python3


import argparse
from dataclasses import dataclass  # builtin in Python 3.7+
import functools
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
  """Overlays efficiencies as a function of one binning variable and stepping through the bins of the other variable given by `steppingVar` for all fits with matching 2D binning"""
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


if __name__ == "__main__":
  makePlots.printGitInfo()
  ROOT.gROOT.SetBatch(True)
  plotTools.setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Overlays efficiency graphs obtained from BruFit results in given directories using given labels.")
  parser.add_argument("--fitResult", type = str, nargs = 2, metavar = ("DIR_PATH", "LABEL"), action = "append", help = "The path to the BruFit output directory of the fit that should be added to the overlay and the corresponding legend label; can be specified multiple times")
  args = parser.parse_args()

  # run = "_041003"
  fitResults: Tuple[Tuple[str, str], ...] = (
    # no extra cuts
    # ("noCut/BruFitOutput.sig_allFixed",        "MC truth"),
    # ("noCut/BruFitOutput.sig_allFudge",        "bggenSig sig fudge"),
    # ("noCut/BruFitOutput.sig_bkg_allFixed",    "bggenSig sig+bkg fixed"),
    # ("noCut/BruFitOutput.sig_bkg_bkgAllFudge", "bggenSig sig+bkg bkg fudge"),
    # ("noCut/BruFitOutput.sig_bkg_allFudge",    "bggenSig sig+bkg all fudge"),
    # ("noCut/BruFitOutput.bggen_allFixed",      "bggen all fixed"),
    # ("noCut/BruFitOutput.bggen_bkgAllFudge",   "bggen sig fixed"),
    # ("noCut/BruFitOutput.bggen_allFudge",      "bggen all free"),
    # ("noCut/BruFitOutput.data_allFixed",       "Data all fixed"),
    # ("noCut/BruFitOutput.data_bkgAllFudge",    "Data sig fixed"),
    # ("noCut/BruFitOutput.data_allFudge",       "Data all free"),
    #
    # no unused showers
    # ("noShowers/BruFitOutput.sig_allFixed",                "MC truth"),
    # # ("noShowers/BruFitOutput.sig_allFudge",                "bggenSig sig fudge"),
    # # ("noShowers/BruFitOutput.sig_bkg_allFixed",            "bggenSig sig+bkg fixed"),
    # # ("noShowers/BruFitOutput.sig_bkg_bkgAllFudge",         "bggenSig sig+bkg bkg fudge"),
    # # ("noShowers/BruFitOutput.sig_bkg_allFudge",            "bggenSig sig+bkg all fudge"),
    # ("noShowers/BruFitOutput.bggen_allFixed",              "bggen fixed"),
    # ("noShowers/BruFitOutput.bggen_bkgAllFudge",           "bggen bkg fugde"),
    # ("noShowers/BruFitOutput.bggen_allFudge",              "bggen all fudge"),
    #
    # ("noShowers/BruFitOutput.sig_allFixed",                      "bggen MC"),
    # (f"noShowers/BruFitOutput.data{run}_allFixed",               "data fixed"),
    # (f"noShowers/BruFitOutput.data{run}_bkgSmear",               "data bkg smear"),
    # # (f"noShowers/BruFitOutput.data{run}_bkgShift",               "data bkg shift"),
    # (f"noShowers/BruFitOutput.data{run}_bkgScale",               "data bkg scale"),
    # # (f"noShowers/BruFitOutput.data{run}_bkgFixSmear",            "data bkg fix smear"),
    # (f"noShowers/BruFitOutput.data{run}_bkgFixShift",            "data bkg fix shift"),
    # # (f"noShowers/BruFitOutput.data{run}_bkgFixScale",            "data bkg fix scale"),
    # (f"noShowers/BruFitOutput.data{run}_bkgAllFudge",            "data bkg fudge"),
    #
    # (f"noShowers/BruFitOutput.data{run}_sigSmear",               "data sig smear"),
    # # (f"noShowers/BruFitOutput.data{run}_sigShift",               "data sig shift"),
    # # (f"noShowers/BruFitOutput.data{run}_sigScale",               "data sig scale"),
    # (f"noShowers/BruFitOutput.data{run}_sigFixSmear",            "data sig fix smear"),
    # # (f"noShowers/BruFitOutput.data{run}_sigFixShift",            "data sig fix shift"),
    # (f"noShowers/BruFitOutput.data{run}_sigFixScale",            "data sig fix scale"),
    # (f"noShowers/BruFitOutput.data{run}_sigShift_bkgShift",      "data sig+bkg shift"),
    # (f"noShowers/BruFitOutput.data{run}_sigShift_bkgShiftScale", "data sig shift bkg shift+scale"),
    # (f"noShowers/BruFitOutput.data{run}_sigAllFudge",            "data sig fudge"),
    #
    # (f"noShowers/BruFitOutput.data{run}_allFudge",               "data all fudge"),
    #
    # (f"noShowers/BruFitOutput.data{run}_sigAllFixed_bkgDoubleGaussian",            "data bkg Double-Gauss"),
    # (f"noShowers/BruFitOutput.data{run}_sigAllFixed_bkgDoubleGaussian_SameMean",   "data bkg Double-Gauss same mean"),
    # (f"noShowers/BruFitOutput.data{run}_sigAllFixed_bkgSkewedGaussian_ExpMod",     "data bkg ExpModGauss"),
    # (f"noShowers/BruFitOutput.data{run}_sigAllFixed_bkgSkewedGaussian_Log",        "data bkg LogNormal"),
    # (f"noShowers/BruFitOutput.data{run}_sigAllFixed_bkgSkewedGaussian_SkewNormal", "data bkg SkewNormal"),
    #
    # (f"noShowers/BruFitOutput.data{run}_sigAllFudge_bkgDoubleGaussian",            "data sig fudge bkg Double-Gauss"),
    # (f"noShowers/BruFitOutput.data{run}_sigAllFudge_bkgDoubleGaussian_SameMean",   "data sig fudge bkg Double-Gauss same mean"),
    # (f"noShowers/BruFitOutput.data{run}_sigAllFudge_bkgSkewedGaussian_ExpMod",     "data sig fudge bkg ExpModGauss"),
    # (f"noShowers/BruFitOutput.data{run}_sigAllFudge_bkgSkewedGaussian_Log",        "data sig fudge bkg LogNormal"),
    # (f"noShowers/BruFitOutput.data{run}_sigAllFudge_bkgSkewedGaussian_SkewNormal", "data sig fudge bkg SkewNormal"),
    #
    # ("2017_01-ver03/noShowers/BruFitOutput.bggen_2017_01-ver03_allFixed", "2017_01-ver03"),
    # ("2018_01-ver02/noShowers/BruFitOutput.bggen_2018_01-ver02_allFixed", "2018_01-ver02"),
    # ("2018_08-ver02/noShowers/BruFitOutput.bggen_2018_08-ver02_allFixed", "2018_08-ver02"),
    # ("2019_11-ver01/noShowers/BruFitOutput.bggen_2019_11-ver01_allFixed", "2019_11-ver01"),
    #
    ("2017_01-ver03/noShowers/BruFitOutput.data_2017_01-ver03_allFixed",  "2017_01-ver03"),
    ("2018_01-ver02/noShowers/BruFitOutput.data_2018_01-ver02_allFixed",  "2018_01-ver02"),
    ("2018_08-ver02/noShowers/BruFitOutput.data_2018_08-ver02_allFixed",  "2018_08-ver02"),
    ("2019_11-ver01/noShowers/BruFitOutput.data_2019_11-ver01_allFixed",  "2019_11-ver01"),
    #
    # ("2017_01-ver03/noShowers/BruFitOutput.bggen_2017_01-ver03_allFixed", "bggen MC"),
    # ("2017_01-ver03/noShowers/BruFitOutput.data_2017_01-ver03_allFixed",  "Real Data"),
    #
    # ("2018_01-ver02/noShowers/BruFitOutput.bggen_2018_01-ver02_allFixed", "bggen MC"),
    # ("2018_01-ver02/noShowers/BruFitOutput.data_2018_01-ver02_allFixed",  "Real Data"),
    #
    # ("2018_08-ver02/noShowers/BruFitOutput.bggen_2018_08-ver02_allFixed", "bggen MC"),
    # ("2018_08-ver02/noShowers/BruFitOutput.data_2018_08-ver02_allFixed",  "Real Data"),
    #
    # ("2019_11-ver01/noShowers/BruFitOutput.bggen_2019_11-ver01_allFixed", "bggen MC"),
    # ("2019_11-ver01/noShowers/BruFitOutput.data_2019_11-ver01_allFixed",  "Real Data"),
  )
  skipBlack = True
  # skipBlack = False
  if args.fitResult:
    fitResultDirNames = tuple(fitResult[0] for fitResult in args.fitResult)
    fitLabels         = tuple(fitResult[1] for fitResult in args.fitResult)
  else:
    # fitResultDirNames = tuple(f"./fits/2018_01-ver02/{fitResult[0]}" for fitResult in fitResults)
    fitResultDirNames = tuple(f"./fits/{fitResult[0]}" for fitResult in fitResults)
    fitLabels         = tuple(fitResult[1] for fitResult in fitResults)
  effInfos, binVarNames = getEfficiencies(fitResultDirNames, fitLabels)
  print("Overlaying efficiencies")
  if effInfos:
    if binVarNames:
      for binningVars in binVarNames:
        if len(binningVars) == 1:
          overlayEfficiencies1D(effInfos, binningVars[0], pdfDirName = "overlays", skipBlack = skipBlack)
        if len(binningVars) == 2:
          overlayEfficiencies2D(effInfos, binningVars[:2], steppingVar = binningVars[1], pdfDirName = "overlays", skipBlack = skipBlack)
