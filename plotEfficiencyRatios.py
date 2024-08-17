#!/usr/bin/env python3


from __future__ import annotations

from collections import defaultdict
from collections.abc import (
  Mapping,
  Sequence,
)
import functools
import os

import ROOT

from overlayEfficiencies import getEfficiencies
from plotBeautifiers import Lines
from plotEfficiencies import (
  EffInfo,
  getEffValuesForGraph1D,
  getEffValuesForGraph2D,
)
from plotFitResults import (
  BinningInfo,
  getAxisInfoForBinningVar,
  plotGraphs1D,
)
import plotTools
from plotTools import (
  calcRatioOfGraphs1D,
  calcRatioOfGraphs2D,
  getGraph1DFromValues,
  getGraph2DFromValues,
  getRangeOfGraph,
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
  effInfos:          Mapping[str, Mapping[tuple[str, str], Sequence[EffInfo]]],  # [ratioLabel][(fitResultDirName, fitLabel)][bin index]
  binningInfo:       BinningInfo,  # 1D binning information
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str        = "Proton_4pi_",
  pdfFileNameSuffix: str        = "",
  graphTitle:        str | None = None,
) -> None:
  """Plots efficiency ratios as a function of `binningVar` for all given fits with 1D binning"""
  assert len(binningInfo.varNames) == 1, f"Need 1-dimensional binning, but got {binningInfo}"
  binningVar = binningInfo.varNames[0]
  print(f"Plotting efficiency ratio for binning variable '{binningVar}'")
  ratioGraphs: list[tuple[str, ROOT.TGraphErrors]] = []
  for ratioLabel, effInfosForLabel in effInfos.items():
    assert len(effInfosForLabel) == 2, f"Expect exactly 2 data samples to calculate ratio; but got {effInfosForLabel}"
    graphs: list[ROOT.TGraphErrors] = []
    for effInfo in effInfosForLabel.values():
      graphs.append(getGraph1DFromValues(getEffValuesForGraph1D(binningVar, effInfo)))
    ratioGraphs.append((ratioLabel, calcRatioOfGraphs1D(graphs)))
  plotGraphs1D(
    graphOrGraphs     = ratioGraphs,
    binningInfo       = binningInfo,
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
  effInfos:          Mapping[str, Mapping[tuple[str, str], Sequence[EffInfo]]],  # [ratioLabel][(fitResultDirName, fitLabel)][bin index]
  binningInfo:       BinningInfo,  # 2D binning information
  steppingVar:       str,  # each 1D graph corresponds to a slice in this variable
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str        = "Proton_4pi_",
  pdfFileNameSuffix: str        = "",
  graphTitle:        str | None = None,
  fitGraphs:         bool       = False,
) -> None:
  """Plots efficiency ratios as a function of one binning variable while stepping through the bins of another variable given by `steppingVar` for all fits with matching 2D binning"""
  assert len(binningInfo.varNames) == 2, f"Need 2-dimensional binning, but got {binningInfo}"
  binningVars = binningInfo.varNames[:2]
  print(f"Plotting efficiency ratios for binning variables '{binningVars}' stepping through bins in '{steppingVar}'")
  assert steppingVar in binningVars, f"Given stepping variable '{steppingVar}' must be in binning variables '{binningVars}'"
  _, steppingVarLabel, steppingVarUnit = getAxisInfoForBinningVar(steppingVar)
  binningVarIndex  = 0 if steppingVar == binningVars[1] else 1
  graphsToOverlay: dict[tuple[float, float], list[tuple[str, ROOT.TGraphErrors]]] = defaultdict(list)
  for ratioLabel, effInfosForLabel in effInfos.items():
    assert len(effInfosForLabel) == 2, f"Expect exactly 2 data samples to calculate ratio; but got {effInfosForLabel}"
    # get efficiencies as 2D graphs
    efficiencyGraphs2D: tuple[ROOT.TGraph2DErrors, ...] = tuple(getGraph2DFromValues(getEffValuesForGraph2D(binningVars, efficiencies))
                                                                for efficiencies in effInfosForLabel.values())
    # calculate 2D graph for efficiency ratios and slice it to 1D graphs
    ratioGraphs1D: dict[tuple[float, float], ROOT.TGraphErrors] = slice2DGraph(
      calcRatioOfGraphs2D(efficiencyGraphs2D, ratioRange = (None, 1.5)),
      plotTools.Graph2DVar.x if steppingVar == binningVars[0] else plotTools.Graph2DVar.y
    )
    for steppingVarBinRange, graph in ratioGraphs1D.items():
      if fitGraphs:
        _, _, xMax, _ = getRangeOfGraph(graph)
        epsilon = 1e-3 * graph.GetErrorX(0)  # used to ensure that right bin is chosen
        graph.Fit("pol0", "SEX0EM", "", 0.8 + epsilon, max(0.8 + epsilon, xMax - epsilon))
      graphsToOverlay[steppingVarBinRange].append((ratioLabel, graph))
  for steppingVarBinRange, graphs in graphsToOverlay.items():
    # overlay 1D graphs for current bin of stepping variable
    steppingVarTitle = f"{steppingVarBinRange[0]} {steppingVarUnit} < {steppingVarLabel} < " \
                       f"{steppingVarBinRange[1]} {steppingVarUnit}"
    plotGraphs1D(
      graphOrGraphs     = graphs,
      binningInfo       = binningInfo.varBinningInfo(binningVars[binningVarIndex]),
      yAxisTitle        = "Efficiency Ratio",
      pdfDirName        = pdfDirName,
      pdfFileBaseName   = "mm2_effratio",
      pdfFileNamePrefix = pdfFileNamePrefix,
      pdfFileNameSuffix = f"_{steppingVar}_{steppingVarBinRange[0]}_{steppingVarBinRange[1]}{pdfFileNameSuffix}",
      graphTitle        = f"{graphTitle}, {steppingVarTitle}",
      graphMinimum      = 0.5,
      graphMaximum      = 1.5,
      skipBlack         = True if len(graphs) > 1 else False,
      drawLegend        = True if len(graphs) > 1 else False,
      # forceXRange       = (1, 6),
      forceXRange       = (0, 9),
      beautifiers       = (HLINES,),
    )


#TODO move to plotEfficiencies and also use in plotEfficiencies2DColzText()
def getHist2DFromEfficiencies(
  binningInfo:   BinningInfo,  # 2D binning information
  efficiencies:  Sequence[EffInfo],
  histName:      str,
) -> ROOT.TH2D:
  """Constructs 2D histogram of efficiency as a function of given binning variables; works only for equidistant binning"""
  assert len(binningInfo.varNames) == 2, f"Need 2-dimensional binning, but got {binningInfo}"
  # TGraph2D always performs interpolation when drawn with COLZ -> construct TH2 with matching binning
  binningVars = tuple(reversed(binningInfo.varNames[:2]))  # swap p and theta axes
  # filter out relevant efficiencies
  effInfos = tuple(
    effInfo for effInfo in efficiencies
    if (binningVars[0] in effInfo.binInfo.varNames) and (binningVars[1] in effInfo.binInfo.varNames) and (len(effInfo.binInfo.varNames) == 2)
  )
  binningVarLabels: list[str] = [""] * len(binningVars)
  binningVarUnits:  list[str] = [""] * len(binningVars)
  for index, binningVar in enumerate(binningVars):
    _, binningVarLabels[index], binningVarUnits[index] = getAxisInfoForBinningVar(binningVar)
  efficiencyHist = ROOT.TH2D(
    histName,
    f";{binningVarLabels[0]} ({binningVarUnits[0]})"
    f";{binningVarLabels[1]} ({binningVarUnits[1]})",
    binningInfo.varNmbBins(binningVars[0]), *binningInfo.varRange(binningVars[0]),
    binningInfo.varNmbBins(binningVars[1]), *binningInfo.varRange(binningVars[1]),
  )
  # fill histogram
  for effInfo in effInfos:
    efficiencyHist.SetBinContent(
      efficiencyHist.FindBin(effInfo.binInfo.center(binningVars[0]), effInfo.binInfo.center(binningVars[1])),
      effInfo.value.nominal_value
    )
  return efficiencyHist


def overlayEfficiencyRatios2DColzText(
  effInfos:          Mapping[str, Mapping[tuple[str, str], Sequence[EffInfo]]],  # [ratioLabel][(fitResultDirName, fitLabel)][bin index]
  binningInfo:       BinningInfo,  # 2D binning information
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNamePrefix: str        = "Proton_4pi_",
  pdfFileNameSuffix: str        = "",
  histTitle:         str | None = None,
) -> None:
  """Plots efficiency ratios as a function of given binning variables for 2-dimensional binning using 'COLZ TEXT' option; works only for equidistant binning"""
  assert len(binningInfo.varNames) == 2, f"Need 2-dimensional binning, but got {binningInfo}"
  binningVars = binningInfo.varNames[:2]
  print(f"Plotting efficiency ratios as a function of binning variables '{binningVars}' assuming equidistant binning")
  for ratioLabel, effInfosForLabel in effInfos.items():
    assert len(effInfosForLabel) == 2, f"Expect exactly 2 data samples to calculate ratio; but got {effInfosForLabel}"
    # get efficiencies as 2D histograms
    effHists2D: tuple[ROOT.TH2D, ...] = tuple(
      getHist2DFromEfficiencies(binningInfo, effs, histName = f"hEfficiency2D_{binningVars[0]}_{binningVars[1]}_{fitLabel}")
      for (_, fitLabel), effs in effInfosForLabel.items()
    )
    # calculate ratio histogram
    ratioHist2D = effHists2D[0].Clone()
    ratioHist2D.Divide(effHists2D[1])
    # draw ratio histogram
    canv = ROOT.TCanvas(f"{pdfFileNamePrefix}mm2_effratio_{binningVars[0]}_{binningVars[1]}_{ratioLabel.replace(' ', '_')}{pdfFileNameSuffix}", "")
    canv.SetRightMargin(0.15)  # make space for color axis
    ROOT.gStyle.SetPalette(ROOT.kVisibleSpectrum)
    if histTitle is not None:
      ratioHist2D.SetTitle(histTitle)
    ratioHist2D.SetZTitle("Efficiency Ratio")
    ratioHist2D.SetMinimum(0.93)  # adjust range such that green part of color scale corresponds to approximately 1.00 +- 0.02
    ratioHist2D.SetMaximum(1.09)
    ROOT.gStyle.SetPaintTextFormat("1.3f")
    ratioHist2D.Draw("COLZ0 TEXT")  # make sure values outside of z range are also plotted
    ratioHist2D.SetStats(False)
    redrawFrame(canv)
    canv.SaveAs(f"{pdfDirName}/{canv.GetName()}_ColzText.pdf")


if __name__ == "__main__":
  printGitInfo()
  ROOT.gROOT.SetBatch(True)
  setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")

  fitRootDir = "./fits.nominal"
  pdfDirName = makeDirPath("./ratios")
  # fitRootDir = "./fits.pionComparison"
  # pdfDirName = makeDirPath("./ratios.pionComparison")
  useMissing = True

  ratiosToPlot: dict[str, tuple[str, str]] = {
    # # all fixed vs. free parameters for data
    # "2017_01-ver03" : (
    #   f"{fitRootDir}/2017_01-ver03/BruFitOutput.data_2017_01-ver03_allFixed",
    #   f"{fitRootDir}/2017_01-ver03/BruFitOutput.data_2017_01-ver03_sigSmear",
    # ),
    # "2018_01-ver02" : (
    #   f"{fitRootDir}/2018_01-ver02/BruFitOutput.data_2018_01-ver02_allFixed",
    #   f"{fitRootDir}/2018_01-ver02/BruFitOutput.data_2018_01-ver02_sigSmear",
    # ),
    # "2018_08-ver02" : (
    #   f"{fitRootDir}/2018_08-ver02/BruFitOutput.data_2018_08-ver02_allFixed",
    #   f"{fitRootDir}/2018_08-ver02/BruFitOutput.data_2018_08-ver02_sigSmear",
    # ),
    # "2019_11-ver01" : (
    #   f"{fitRootDir}/2019_11-ver01/BruFitOutput.data_2019_11-ver01_allFixed",
    #   f"{fitRootDir}/2019_11-ver01/BruFitOutput.data_2019_11-ver01_sigSmear",
    # ),
    # old data vs. MC
    "2017_01-ver03" : (
      f"{fitRootDir}/2017_01-ver03/BruFitOutput.data_2017_01-ver03_allFixed",
      f"{fitRootDir}/2017_01-ver03/BruFitOutput.bggen_2017_01-ver03_allFixed",
    ),
    # "2018_01-ver02" : (
    #   f"{fitRootDir}/2018_01-ver02/BruFitOutput.data_2018_01-ver02_allFixed",
    #   f"{fitRootDir}/2018_01-ver02/BruFitOutput.bggen_2018_01-ver02_allFixed",
    # ),
    # "2018_08-ver02" : (
    #   f"{fitRootDir}/2018_08-ver02/BruFitOutput.data_2018_08-ver02_allFixed",
    #   f"{fitRootDir}/2018_08-ver02/BruFitOutput.bggen_2018_08-ver02_allFixed",
    # ),
    # "2019_11-ver01" : (
    #   f"{fitRootDir}/2019_11-ver01/BruFitOutput.data_2019_11-ver01_allFixed",
    #   f"{fitRootDir}/2019_11-ver01/BruFitOutput.bggen_2019_11-ver01_allFixed",
    # ),
    # # comparison old vs. good ToF
    # "2018_01 MC" : (
    #   f"{fitRootDir}/2018_01-ver02_goodToF/BruFitOutput.bggen_2018_01-ver02_goodToF_allFixed",
    #   f"{fitRootDir}/2018_01-ver02/BruFitOutput.bggen_2018_01-ver02_allFixed",
    # ),
    # "2018_01 Data" : (
    #   f"{fitRootDir}/2018_01-ver02_goodToF/BruFitOutput.data_2018_01-ver02_goodToF_allFixed",
    #   f"{fitRootDir}/2018_01-ver02/BruFitOutput.data_2018_01-ver02_allFixed",
    # ),
    # "2018_08 MC" : (
    #   f"{fitRootDir}/2018_08-ver02_goodToF/BruFitOutput.bggen_2018_08-ver02_goodToF_allFixed",
    #   f"{fitRootDir}/2018_08-ver02/BruFitOutput.bggen_2018_08-ver02_allFixed",
    # ),
    # "2018_08 Data" : (
    #   f"{fitRootDir}/2018_08-ver02_goodToF/BruFitOutput.data_2018_08-ver02_goodToF_allFixed",
    #   f"{fitRootDir}/2018_08-ver02/BruFitOutput.data_2018_08-ver02_allFixed",
    # ),
    # # good ToF data vs. MC
    # "2018_01 good ToF" : (
    #   f"{fitRootDir}/2018_01-ver02_goodToF/BruFitOutput.data_2018_01-ver02_goodToF_allFixed",
    #   f"{fitRootDir}/2018_01-ver02_goodToF/BruFitOutput.bggen_2018_01-ver02_goodToF_allFixed",
    # ),
    # "2018_08 good ToF" : (
    #   f"{fitRootDir}/2018_08-ver02_goodToF/BruFitOutput.data_2018_08-ver02_goodToF_allFixed",
    #   f"{fitRootDir}/2018_08-ver02_goodToF/BruFitOutput.bggen_2018_08-ver02_goodToF_allFixed",
    # ),
  }

  # title = "Real Data all fixed / sig smear"
  title = "Real Data / MC"
  # title = "Good ToF / All Files"
  effInfos:     dict[str, dict[tuple[str, str], list[EffInfo]]] = {}
  binningInfos: dict[str, list[BinningInfo | None]]             = {}
  for ratioLabel, fitResults in ratiosToPlot.items():
    effInfos[ratioLabel], binningInfos[ratioLabel] = getEfficiencies(
      fitResultDirNames = tuple(fitResult for fitResult in fitResults),
      useMissing = useMissing
    )
  print("Plotting efficiency ratios")
  if effInfos and binningInfos:
    firstBinningInfos = next(iter(binningInfos.values()))  # get first entry in binningInfos
    for ratioLabel, binningInfosForLabel in binningInfos.items():
      assert binningInfosForLabel == firstBinningInfos, f"Data samples have different binnings: '{ratioLabel}' = {binningInfosForLabel} vs. '{next(iter(binningInfos.keys()))}' = {firstBinningInfos}"
    if firstBinningInfos:
      for binningInfo in firstBinningInfos:
        if binningInfo:
          if len(binningInfo.varNames) == 1:
            overlayEfficiencyRatios1D(effInfos, binningInfo, pdfDirName, graphTitle = title)
          elif len(binningInfo.varNames) == 2:
            overlayEfficiencyRatios2DSlices  (effInfos, binningInfo, steppingVar = binningInfo.varNames[1], pdfDirName = pdfDirName, graphTitle = title)
            overlayEfficiencyRatios2DColzText(effInfos, binningInfo, pdfDirName = pdfDirName, histTitle = title)
