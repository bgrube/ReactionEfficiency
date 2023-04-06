#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


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
  pdfFileNameSuffix: str   = "",
  particle:          str   = "Proton",
  channel:           str   = "4pi",
  markerSize:        float = 0.75,
):
  '''Overlays efficiencies as a function of `binningVar` for all given fits'''
  print(f"Overlaying efficiencies for binning variable '{binningVar}'")
  efficiencyMultiGraph = ROOT.TMultiGraph()  # type: ignore
  efficiencyGraphs = {}  # store graphs here to keep them in memory
  shiftFraction = 0
  color = 1
  for (fitResultDirName, fitLabel), efficiencies in effInfos.items():
    graph = efficiencyGraphs[fitResultDirName] = plotFitResults.getParValueGraph1D(plotEfficiencies.getEffValuesForGraph1D(binningVar, efficiencies), shiftFraction)
    graph.SetTitle(fitLabel)            # type: ignore
    graph.SetMarkerStyle(ROOT.kCircle)  # type: ignore
    graph.SetMarkerSize(markerSize)     # type: ignore
    graph.SetMarkerColor(color)         # type: ignore
    graph.SetLineColor(color)           # type: ignore
    shiftFraction += 0.01
    color += 1
    efficiencyMultiGraph.Add(graph)
  efficiencyMultiGraph.SetTitle(f"{particle} Track-Finding Efficiency ({channel})")
  assert binningVar in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVar}'"
  efficiencyMultiGraph.GetXaxis().SetTitle(f"{BINNING_VAR_PLOT_INFO[binningVar]['label']} ({BINNING_VAR_PLOT_INFO[binningVar]['unit']})")
  efficiencyMultiGraph.GetYaxis().SetTitle("Efficiency")
  efficiencyMultiGraph.SetMinimum(0)
  efficiencyMultiGraph.SetMaximum(1)
  fitLabels = tuple(key[1] for key in effInfos.keys())
  canv = ROOT.TCanvas(f"{particle}_{channel}_mm2_eff_{binningVar}_{'_'.join(fitLabels)}{pdfFileNameSuffix}", "")  # type: ignore
  efficiencyMultiGraph.Draw("APZ")
  canv.BuildLegend()
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)  # type: ignore
  makePlots.setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")  # type: ignore

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Overlays efficiency graphs obtained from BruFit results in given directories.")
  parser.add_argument("fitResultDirNames", type = str, nargs = "+", default = ["./BruFitOutput"], help = "The paths to the BruFit output directories that should be overlaid; (default: %(default)s)")
  parser.add_argument("--fitLabels", type = str, nargs = "*", default = [], help = "The legend labels for each fit (same order as directory names); (default: directory names)")
  args = parser.parse_args()
  # fitVariable = "MissingMassSquared_Measured"

  fitResultDirNames = ["./BruFitOutput.data_allFixed", "./BruFitOutput.data_bkgAllFudge", "./BruFitOutput.data_allFudge"]
  fitLabels         = ["All fixed", "Sig fixed", "All free"]
  effInfos, binVarNames = getEfficiencies(fitResultDirNames, fitLabels)
  print("Overlaying efficiencies")
  print(effInfos)
  if effInfos:
    if binVarNames:
      for binningVars in binVarNames:
        if len(binningVars) == 1:
          overlayEfficiencies(effInfos, binningVars[0], "foo")
