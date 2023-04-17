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
  assert len(fitResultDirNames) == len(set(fitResultDirNames)), f"list of fit-result directory names '{fitResultDirNames}' must have unique elements"
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

  # # echo and parse command line
  # print(f"Script was called using: '{' '.join(sys.argv)}'")
  # parser = argparse.ArgumentParser(description="Overlays efficiency graphs obtained from BruFit results in given directories.")
  # parser.add_argument("fitResultDirNames", type = str, nargs = "+", default = ["./BruFitOutput"], help = "The paths to the BruFit output directories that should be overlaid; (default: %(default)s)")
  # parser.add_argument("--fitLabels", type = str, nargs = "*", default = [], help = "The legend labels for each fit (same order as directory names); (default: directory names)")
  # args = parser.parse_args()

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
    # no unused showers
    ("noShowers/BruFitOutput.sig_allFixed",                "NUS MC truth"),
    # ("noShowers/BruFitOutput.sig_allFudge",                "NUS bggenSig sig fudge"),
    # ("noShowers/BruFitOutput.sig_bkg_allFixed",            "NUS bggenSig sig+bkg fixed"),
    # ("noShowers/BruFitOutput.sig_bkg_bkgAllFudge",         "NUS bggenSig sig+bkg bkg fudge"),
    # ("noShowers/BruFitOutput.sig_bkg_allFudge",            "NUS bggenSig sig+bkg all fudge"),
    # ("noShowers/BruFitOutput.bggen_allFixed",              "NUS bggen fixed"),
    # ("noShowers/BruFitOutput.bggen_bkgAllFudge",           "NUS bggen bkg fugde"),
    # ("noShowers/BruFitOutput.bggen_allFudge",              "NUS bggen all fudge"),
    #
    ("noShowers/BruFitOutput.data_allFixed",               "NUS data fixed"),
    # ("noShowers/BruFitOutput.data_bkgSmear",               "NUS data bkg smear"),
    # # ("noShowers/BruFitOutput.data_bkgShift",               "NUS data bkg shift"),
    # ("noShowers/BruFitOutput.data_bkgScale",               "NUS data bkg scale"),
    # # ("noShowers/BruFitOutput.data_bkgFixSmear",            "NUS data bkg fix smear"),
    # ("noShowers/BruFitOutput.data_bkgFixShift",            "NUS data bkg fix shift"),
    # # ("noShowers/BruFitOutput.data_bkgFixScale",            "NUS data bkg fix scale"),
    # ("noShowers/BruFitOutput.data_bkgAllFudge",            "NUS data bkg fudge"),
    #
    # ("noShowers/BruFitOutput.data_sigSmear",               "NUS data sig smear"),
    # # ("noShowers/BruFitOutput.data_sigShift",               "NUS data sig shift"),
    # # ("noShowers/BruFitOutput.data_sigScale",               "NUS data sig scale"),
    # ("noShowers/BruFitOutput.data_sigFixSmear",            "NUS data sig fix smear"),
    # # ("noShowers/BruFitOutput.data_sigFixShift",            "NUS data sig fix shift"),
    # ("noShowers/BruFitOutput.data_sigFixScale",            "NUS data sig fix scale"),
    # ("noShowers/BruFitOutput.data_sigShift_bkgShift",      "NUS data sig+bkg shift"),
    # ("noShowers/BruFitOutput.data_sigShift_bkgShiftScale", "NUS data sig shift bkg shift+scale"),
    # ("noShowers/BruFitOutput.data_allFudge",               "NUS data all fudge"),
    #
    # ("noShowers/BruFitOutput.data_sigAllFixed_bkgDoubleGaussian",            "NUS data bkg Double-Gauss"),
    # ("noShowers/BruFitOutput.data_sigAllFixed_bkgDoubleGaussian_SameMean",   "NUS data bkg Double-Gauss same mean"),
    # ("noShowers/BruFitOutput.data_sigAllFixed_bkgSkewedGaussian_ExpMod",     "NUS data bkg ExpModGauss"),
    # ("noShowers/BruFitOutput.data_sigAllFixed_bkgSkewedGaussian_Log",        "NUS data bkg LogNormal"),
    # ("noShowers/BruFitOutput.data_sigAllFixed_bkgSkewedGaussian_SkewNormal", "NUS data bkg SkewNormal"),
    #
    ("noShowers/BruFitOutput.data_sigAllFudge_bkgDoubleGaussian",            "NUS data sig fudge bkg Double-Gauss"),
    ("noShowers/BruFitOutput.data_sigAllFudge_bkgDoubleGaussian_SameMean",   "NUS data sig fudge bkg Double-Gauss same mean"),
    ("noShowers/BruFitOutput.data_sigAllFudge_bkgSkewedGaussian_ExpMod",     "NUS data sig fudge bkg ExpModGauss"),
    ("noShowers/BruFitOutput.data_sigAllFudge_bkgSkewedGaussian_Log",        "NUS data sig fudge bkg LogNormal"),
    ("noShowers/BruFitOutput.data_sigAllFudge_bkgSkewedGaussian_SkewNormal", "NUS data sig fudge bkg SkewNormal"),
  ]
  fitResultDirNames = tuple(f"./fits/{fitResult[0]}" for fitResult in fitResults)
  fitLabels         = tuple(fitResult[1] for fitResult in fitResults)
  effInfos, binVarNames = getEfficiencies(fitResultDirNames, fitLabels)
  print("Overlaying efficiencies")
  if effInfos:
    if binVarNames:
      for binningVars in binVarNames:
        if len(binningVars) == 1:
          overlayEfficiencies(effInfos, binningVars[0], "overlays")
