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
from plotBeautifiers import Lines
from plotEfficiencies import (
  BinInfo,
  EffInfo,
  getEffValuesForGraph1D,
  getEffValuesForGraph2D,
)
from plotFitResults import BINNING_VAR_PLOT_INFO, plotGraphs1D
from plotTools import (
  calcRatioOfGraphs1D,
  calcRatioOfGraphs2D,
  getGraph1DFromValues,
  getGraph2DFromValues,
  getRangeOfGraph,
  Graph2DVar,
  makeDirPath,
  printGitInfo,
  redrawFrame,
  setupPlotStyle,
  slice2DGraph,
)


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


# horizontal lines in graphs
HLINES = Lines(defaultColor = ROOT.kGray + 1, orientation = Lines.Orientation.horizontal, drawContentOverLines = True).set((0.95, 1.00, 1.05))


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
      graphs.append(getGraph1DFromValues(getEffValuesForGraph1D(binningVar, effInfo)))
    ratioGraphs.append((ratioLabel, calcRatioOfGraphs1D(graphs)))
  plotGraphs1D(
    graphOrGraphs     = ratioGraphs,
    binningVar        = binningVar,
    yAxisTitle        = "Efficiency Ratio",
    pdfDirName        = pdfDirName,
    pdfFileBaseName   = "mm2_effratio",
    pdfFileNamePrefix = pdfFileNamePrefix,
    pdfFileNameSuffix = pdfFileNameSuffix,
    graphTitle        = f"Efficiency Ratio" if graphTitle is None else graphTitle,
    graphMinimum      = 0.5,
    graphMaximum      = 1.5,
    skipBlack         = True if len(ratioGraphs) > 1 else False,
    drawLegend        = True if len(ratioGraphs) > 1 else False,
    beautifiers       = (HLINES,),
  )


def overlayEfficiencyRatios2DSlices(
  effInfos:          Mapping[str, Mapping[Tuple[str, str], Sequence[EffInfo]]],  # [ratioLabel][(fitResultDirName, fitLabel)][bin index]
  binningVars:       Sequence[str],
  steppingVar:       str,
  pdfDirName:        str,
  pdfFileNamePrefix: str = "Proton_4pi_",
  pdfFileNameSuffix: str = "",
  graphTitle:        Optional[str] = None,
  fitGraphs:         bool = False,
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
    efficiencyGraphs2D: Tuple[ROOT.TGraph2DErrors, ...] = tuple(getGraph2DFromValues(getEffValuesForGraph2D(binningVars, efficiencies))
                                                                for efficiencies in effInfosForLabel.values())
    # calculate 2D graph for efficiency ratios and slice it to 1D graphs
    ratioGraphs1D: Dict[Tuple[float, float], ROOT.TGraphErrors] = slice2DGraph(
      calcRatioOfGraphs2D(efficiencyGraphs2D, ratioRange = (None, 1.5)),
      Graph2DVar.x if steppingVar == binningVars[0] else Graph2DVar.y
    )
    for steppingVarBinRange, graph in ratioGraphs1D.items():
      if fitGraphs:
        _, _, xMax, _ = getRangeOfGraph(graph)
        epsilon = 1e-3 * graph.GetErrorX(0)  # used to ensure that right bin is chosen
        graph.Fit("pol0", "SEX0EM", "", 0.8 + epsilon, max(0.8 + epsilon, xMax - epsilon))
      graphsToOverlay[steppingVarBinRange].append((ratioLabel, graph))
  for steppingVarBinRange, graphs in graphsToOverlay.items():
    # overlay 1D graphs for current bin of stepping variable
    steppingVarLabel = f"{steppingVarBinRange[0]} {BINNING_VAR_PLOT_INFO[steppingVar]['unit']} " \
      f"< {BINNING_VAR_PLOT_INFO[steppingVar]['label']} " \
      f"< {steppingVarBinRange[1]} {BINNING_VAR_PLOT_INFO[steppingVar]['unit']}"
    plotGraphs1D(
      graphs,
      binningVars[binningVarIndex],
      yAxisTitle        = "Efficiency Ratio",
      pdfDirName        = pdfDirName,
      pdfFileBaseName   = "mm2_effratio",
      pdfFileNamePrefix = pdfFileNamePrefix,
      pdfFileNameSuffix = f"_{steppingVar}_{steppingVarBinRange[0]}_{steppingVarBinRange[1]}{pdfFileNameSuffix}",
      graphTitle        = f"{graphTitle}, {steppingVarLabel}",
      graphMinimum      = 0.5,
      graphMaximum      = 1.5,
      skipBlack         = True if len(graphs) > 1 else False,
      drawLegend        = True if len(graphs) > 1 else False,
      forceXRange       = (1, 6),
      beautifiers       = (HLINES,),
    )


#TODO move to plotEfficiencies and also use in plotEfficiencies2DColzText()
def getHist2DFromEfficiencies(
  binningVarsIn: Sequence[str],  # names of binning variables to plot
  efficiencies:  Sequence[EffInfo],
  histName:      str,
) -> ROOT.TH2D:
  """Constructs 2D histogram of efficiency as a function of given binning variables; works only for equidistant binning"""
  # TGraph2D always performs interpolation when drawn with COLZ -> construct TH2 with matching binning
  binningVars = tuple(reversed(binningVarsIn))  # swap p and theta axes
  # filter out relevant efficiencies
  effInfos = tuple(effInfo for effInfo in efficiencies
                   if (binningVars[0] in effInfo.binInfo.varNames) and (binningVars[1] in effInfo.binInfo.varNames) and (len(effInfo.binInfo.varNames) == 2))
  # determine equidistant binning and create histogram
  xWidths = set(effInfo.binInfo.widths[binningVars[0]] for effInfo in effInfos)
  yWidths = set(effInfo.binInfo.widths[binningVars[1]] for effInfo in effInfos)
  assert (len(xWidths) == 1) and (len(yWidths) == 1), f"Binning is not equidistant: x bin widths = {xWidths}; y bin widths = {yWidths}"
  xWidth = tuple(xWidths)[0]
  yWidth = tuple(yWidths)[0]
  xCenters = sorted(set(effInfo.binInfo.centers[binningVars[0]] for effInfo in effInfos))
  yCenters = sorted(set(effInfo.binInfo.centers[binningVars[1]] for effInfo in effInfos))
  xRange = (xCenters[0] - xWidth / 2.0, xCenters[-1] + xWidth / 2.0)
  yRange = (yCenters[0] - yWidth / 2.0, yCenters[-1] + yWidth / 2.0)
  assert binningVars[0] in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVars[0]}'"
  assert binningVars[1] in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVars[1]}'"
  efficiencyHist = ROOT.TH2D(
    histName, f";{BINNING_VAR_PLOT_INFO[binningVars[0]]['label']} ({BINNING_VAR_PLOT_INFO[binningVars[0]]['unit']})"
              f";{BINNING_VAR_PLOT_INFO[binningVars[1]]['label']} ({BINNING_VAR_PLOT_INFO[binningVars[1]]['unit']})",
    len(xCenters), *xRange, len(yCenters), *yRange)
  # fill histogram
  for effInfo in effInfos:
    efficiencyHist.SetBinContent(efficiencyHist.FindBin(effInfo.binInfo.centers[binningVars[0]], effInfo.binInfo.centers[binningVars[1]]),
                                 effInfo.value.nominal_value)
  return efficiencyHist


def overlayEfficiencyRatios2DColzText(
  effInfos:          Mapping[str, Mapping[Tuple[str, str], Sequence[EffInfo]]],  # [ratioLabel][(fitResultDirName, fitLabel)][bin index]
  binningVars:       Sequence[str],
  pdfDirName:        str,
  pdfFileNamePrefix: str = "Proton_4pi_",
  pdfFileNameSuffix: str = "",
  histTitle:         Optional[str] = None,
) -> None:
  """Plots efficiency ratios as a function of given binning variables for 2-dimensional binning using 'COLZ TEXT' option; works only for equidistant binning"""
  print(f"Plotting efficiency ratios as a function of binning variables '{binningVars}' assuming equidistant binning")
  for ratioLabel, effInfosForLabel in effInfos.items():
    assert len(effInfosForLabel) == 2, f"Expect exactly 2 data samples to calculate ratio; but got {effInfosForLabel}"
    # get efficiencies as 2D histograms
    effHists2D: Tuple[ROOT.TH2D, ...] = tuple(getHist2DFromEfficiencies(binningVars, effs, histName = f"hEfficiency2D_{binningVars[0]}_{binningVars[1]}_fitLabel")
                                              for (_, fitLabel), effs in effInfosForLabel.items())
    # calculate ratio histogram
    ratioHist2D = effHists2D[0].Clone()
    ratioHist2D.Divide(effHists2D[1])
    # draw ratio histogram
    canv = ROOT.TCanvas(f"{pdfFileNamePrefix}mm2_effratio_{binningVars[0]}_{binningVars[1]}_{ratioLabel}{pdfFileNameSuffix}", "")
    canv.SetRightMargin(0.15)
    ROOT.gStyle.SetPalette(ROOT.kVisibleSpectrum)
    if histTitle is not None:
      ratioHist2D.SetTitle(histTitle)
    ratioHist2D.SetZTitle("Efficiency Ratio")
    ratioHist2D.SetMinimum(0.93)
    ratioHist2D.SetMaximum(1.09)
    ROOT.gStyle.SetPaintTextFormat("1.3f")
    ratioHist2D.Draw("COLZ TEXT")
    ratioHist2D.SetStats(False)
    redrawFrame(canv)
    canv.SaveAs(f"{pdfDirName}/{canv.GetName()}_ColzText.pdf")


if __name__ == "__main__":
  printGitInfo()
  ROOT.gROOT.SetBatch(True)
  setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")

  # fitRootDir = "./fits"
  # pdfDirName = makeDirPath("./ratios")
  # fitRootDir = "./fits.pionComparison"
  # pdfDirName = makeDirPath("./ratios.pionComparison")
  fitRootDir = "./fits.pionComparison2"
  pdfDirName = makeDirPath("./ratios.pionComparison2")

  ratiosToPlot: Dict[str, Tuple[str, str]] = {
    # "2017_01-ver03" : (
    #   f"{fitRootDir}/2017_01-ver03/noShowers/BruFitOutput.data_2017_01-ver03_allFixed",
    #   f"{fitRootDir}/2017_01-ver03/noShowers/BruFitOutput.bggen_2017_01-ver03_allFixed",
    # ),
    # "2018_01-ver02" : (
    #   f"{fitRootDir}/2018_01-ver02/noShowers/BruFitOutput.data_2018_01-ver02_allFixed",
    #   f"{fitRootDir}/2018_01-ver02/noShowers/BruFitOutput.bggen_2018_01-ver02_allFixed",
    # ),
    # "2018_08-ver02" : (
    #   f"{fitRootDir}/2018_08-ver02/noShowers/BruFitOutput.data_2018_08-ver02_allFixed",
    #   f"{fitRootDir}/2018_08-ver02/noShowers/BruFitOutput.bggen_2018_08-ver02_allFixed",
    # ),
    # "2018_08 good ToF" : (
    #   f"{fitRootDir}/2018_08-ver02_goodToF/noShowers/BruFitOutput.data_2018_08-ver02_goodToF_allFixed",
    #   f"{fitRootDir}/2018_08-ver02_goodToF/noShowers/BruFitOutput.bggen_2018_08-ver02_goodToF_allFixed",
    # ),
    # "2019_11-ver01" : (
    #   f"{fitRootDir}/2019_11-ver01/noShowers/BruFitOutput.data_2019_11-ver01_allFixed",
    #   f"{fitRootDir}/2019_11-ver01/noShowers/BruFitOutput.bggen_2019_11-ver01_allFixed",
    # ),
    "2018_08 MC" : (
      f"{fitRootDir}/2018_08-ver02_goodToF/noShowers/BruFitOutput.bggen_2018_08-ver02_goodToF_allFixed",
      f"{fitRootDir}/2018_08-ver02/noShowers/BruFitOutput.bggen_2018_08-ver02_allFixed",
    ),
    "2018_08 Data" : (
      f"{fitRootDir}/2018_08-ver02_goodToF/noShowers/BruFitOutput.data_2018_08-ver02_goodToF_allFixed",
      f"{fitRootDir}/2018_08-ver02/noShowers/BruFitOutput.data_2018_08-ver02_allFixed",
    ),
  }

  # title = "Real Data / MC"
  title = "Good ToF / All Files"
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
          overlayEfficiencyRatios1D(effInfos, binningVars[0], pdfDirName, graphTitle = title)
        if len(binningVars) == 2:
          overlayEfficiencyRatios2DSlices  (effInfos, binningVars[:2], steppingVar = binningVars[1], pdfDirName = pdfDirName, graphTitle = title)
          overlayEfficiencyRatios2DColzText(effInfos, binningVars[:2], pdfDirName = pdfDirName, histTitle = title)
