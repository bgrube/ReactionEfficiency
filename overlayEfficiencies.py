#!/usr/bin/env python3


import argparse
from dataclasses import dataclass  # builtin in Python 3.7+
import functools
import os
import sys
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple, Union

from uncertainties import UFloat, ufloat

import ROOT

import makePlots
import plotFitResults
from plotFitResults import ParInfo, BINNING_VAR_PLOT_INFO
import plotEfficiencies
from plotEfficiencies import EffInfo


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


def getEfficiencies(
  fitResultDirNames: Sequence[str],
  fitLabels:         Sequence[str] = [],
  dataSets:          Iterable[str] = ["Total", "Found", "Missing"],
) -> Tuple[Dict[Tuple[str, str], List[EffInfo]], Union[List[Tuple[str, ...]], None]]:
  '''Reads yields from given fit directories and calculates efficiencies for each directory'''
  assert len(fitResultDirNames) == len(set(fitResultDirNames)), f"list of fit-result directory names '{fitResultDirNames}' must consist of unique elements"
  assert (not fitLabels) or len(fitLabels) == len(fitResultDirNames), f"Number of given fit labels ({len(fitLabels)}) does not match number of fit directories ({len(fitResultDirNames)})"
  print("Reading yields and calculating efficiencies")
  effInfos:    Dict[Tuple[str, str], List[EffInfo]] = {}    # effInfos[(<fit directory>, <fit label>)][<bin>]
  binVarNames: Union[List[Tuple[str, ...]], None]   = None  # binning variables for each binning
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


def overlayEfficiencies(
  effInfos:          Mapping[Tuple[str, str], Sequence[EffInfo]],
  binningVar:        str,
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNameSuffix: str = "",
  particle:          str = "Proton",
  channel:           str = "4pi",
):
  '''Overlays efficiencies as a function of `binningVar` for all given fits'''
  print(f"Overlaying efficiencies for binning variable '{binningVar}'")
  efficiencyMultiGraph = ROOT.TMultiGraph()  # type: ignore
  efficiencyGraphs = {}  # store graphs here to keep them in memory
  shiftFraction = 0
  styleIndex = 0
  for (fitResultDirName, fitLabel), efficiencies in effInfos.items():
    graph = efficiencyGraphs[fitResultDirName] = plotFitResults.getParValueGraph1D(plotEfficiencies.getEffValuesForGraph1D(binningVar, efficiencies), shiftFraction)
    shiftFraction += 0.01
    graph.SetTitle(fitLabel)  # type: ignore
    makePlots.setCbFriendlyStyle(graph, styleIndex, skipBlack = False if len(effInfos) == 1 else True)
    styleIndex += 1
    efficiencyMultiGraph.Add(graph)
  efficiencyMultiGraph.SetTitle(f"{particle} Track-Finding Efficiency ({channel})")
  assert binningVar in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVar}'"
  efficiencyMultiGraph.GetXaxis().SetTitle(f"{BINNING_VAR_PLOT_INFO[binningVar]['label']} ({BINNING_VAR_PLOT_INFO[binningVar]['unit']})")
  efficiencyMultiGraph.GetYaxis().SetTitle("Efficiency")
  efficiencyMultiGraph.SetMinimum(0)
  efficiencyMultiGraph.SetMaximum(1)
  fitLabels = tuple(key[1].replace(' ', '_') for key in effInfos.keys())
  canv = ROOT.TCanvas(f"{particle}_{channel}_mm2_eff_{binningVar}_{'_'.join(fitLabels)}{pdfFileNameSuffix}", "")  # type: ignore
  efficiencyMultiGraph.Draw("APZ")
  canv.BuildLegend()
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


if __name__ == "__main__":
  makePlots.printGitInfo()
  ROOT.gROOT.SetBatch(True)  # type: ignore
  makePlots.setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")  # type: ignore

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Overlays efficiency graphs obtained from BruFit results in given directories using given labels.")
  parser.add_argument("--fitResult", type = str, nargs = 2, metavar = ("DIR_PATH", "LABEL"), action = "append", help = "The path to the BruFit output directory of the fit that should be added to the overlay and the corresponding legend label; can be specified multiple times")
  args = parser.parse_args()

  # run = "_041003"
  fitResults = [
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
    # # 2018 runs
    # (f"noShowers/BruFitOutput.data_041003_allFixed", "Run 41003"),
    # (f"noShowers/BruFitOutput.data_042030_allFixed", "Run 42030"),
    # (f"noShowers/BruFitOutput.data_042550_allFixed", "Run 42550"),
    # 2020 runs
    (f"noShowers/BruFitOutput.data_071592_allFixed", "Run 71592"),
    (f"noShowers/BruFitOutput.data_071593_allFixed", "Run 71593"),
    (f"noShowers/BruFitOutput.data_071594_allFixed", "Run 71594"),
    (f"noShowers/BruFitOutput.data_071596_allFixed", "Run 71596"),
  ]
  if args.fitResult:
    fitResultDirNames = [fitResult[0] for fitResult in args.fitResult]
    fitLabels         = [fitResult[1] for fitResult in args.fitResult]
  else:
    # fitResultDirNames = tuple(f"./fits/2018_01-ver02/{fitResult[0]}" for fitResult in fitResults)
    fitResultDirNames = tuple(f"./fits/2019_11-ver01/{fitResult[0]}" for fitResult in fitResults)
    fitLabels         = tuple(fitResult[1] for fitResult in fitResults)
  effInfos, binVarNames = getEfficiencies(fitResultDirNames, fitLabels)
  print("Overlaying efficiencies")
  if effInfos:
    if binVarNames:
      for binningVars in binVarNames:
        if len(binningVars) == 1:
          overlayEfficiencies(effInfos, binningVars[0], "overlays")
