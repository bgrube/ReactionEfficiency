#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import argparse
from dataclasses import dataclass  # builtin in Python 3.7+
import functools
import math
import numpy as np
import os
import sys
from typing import Dict, List, Tuple, Union

from uncertainties import UFloat, ufloat

import ROOT

import makePlots
import plotFitResults
from plotFitResults import BinInfo, BinningInfo, ParInfo, BINNING_VAR_PLOT_INFO


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


YIELD_PAR_NAMES = {
  "Signal"     : "Yld_SigPdf",
  "Background" : "Yld_BkgPdf",
}


def readDataIntegralFromFitFile(
  binInfo:     BinInfo,
  fitVariable: str,
) -> ParInfo:
  '''Reads data histogram for given kinematic bin and returns its integral assuming Poissonian uncertainty'''
  fitResultFileName = binInfo.fitResultFileName
  print(f"Reading fitted data histogram from file '{fitResultFileName}'")
  fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")  # type: ignore
  canvName = f"{'' if binInfo.name == 'Overall' else binInfo.name}_{fitVariable}"
  canv = fitResultFile.Get(canvName)
  # get data histogram
  dataFitPad = canv.GetListOfPrimitives().FindObject(f"{canvName}_1")
  dataHist = dataFitPad.GetListOfPrimitives().FindObject("h_DataEvents")
  #TODO reimplement in numpy
  yVals = dataHist.GetY()
  integral = 0
  for i in range(dataHist.GetN()):
    integral += yVals[i]
  fitResultFile.Close()
  return ParInfo(binInfo, {"Signal" : ufloat(integral, math.sqrt(integral))})


def readYieldInfosForBinning(
  binningInfo:   BinningInfo,
  fitVariable:   str,
  readIntegrals: bool = False,
) -> List[ParInfo]:
  '''Reads yields (or histogram integrals if readIntegrals is set) from fit-result files for given binning'''
  yieldInfos = []
  # read overall yields
  overallBinInfo = BinInfo("Overall", {}, binningInfo.dirName)  # special bin for overall fit results
  if os.path.isfile(overallBinInfo.fitResultFileName):
    yieldInfo = readDataIntegralFromFitFile(overallBinInfo, fitVariable) if readIntegrals \
      else plotFitResults.readParInfoForBin(overallBinInfo, YIELD_PAR_NAMES)
    print(f"Read overall yields: {yieldInfo}")
    yieldInfos.append(yieldInfo)  # first entry always contains overall yields
  # read yields for each bin
  for binInfo in binningInfo.infos:
    yieldInfo = readDataIntegralFromFitFile(binInfo, fitVariable) if readIntegrals \
      else plotFitResults.readParInfoForBin(binInfo, YIELD_PAR_NAMES)
    print(f"Read yields for kinematic bin: {yieldInfo}")
    yieldInfos.append(yieldInfo)
  return yieldInfos


@dataclass
class EffInfo:
  '''Stores information about efficiency in a single kinematic bin'''
  binInfo: BinInfo  # info for the kinematic bin
  value:   UFloat   # efficiency value


def calculateEfficiencies(yieldInfos: Dict[str, List[ParInfo]]) -> List[EffInfo]:
  '''Calculates efficiencies from yields'''
  assert len(yieldInfos["Total"]) == len(yieldInfos["Found"]) == len(yieldInfos["Missing"]), "Datasets have different number of kinematic bins"
  effInfos = []
  for index, yieldInfoFound in enumerate(yieldInfos["Found"]):
    yieldInfoMissing = yieldInfos["Missing"][index]
    assert yieldInfoFound.binInfo.isSameBinAs(yieldInfoMissing.binInfo), f"Bin infos for 'Found' and 'Missing' are not identical: {yieldInfoFound.binInfo} vs. {yieldInfoMissing.binInfo}"
    # calculate efficiency
    effInfo = EffInfo(yieldInfoFound.binInfo,
      yieldInfoFound.values["Signal"] / (yieldInfoFound.values["Signal"] + yieldInfoMissing.values["Signal"]))  # type: ignore
    print(f"Effiency = {effInfo}")
    effInfos.append(effInfo)
  return effInfos


def getEffValuesForGraph1D(
  binVarName: str,  # name of the binning variable, i.e. x-axis
  effInfos:   List[EffInfo],
) -> List[Tuple[float, UFloat]]:
  '''Extracts information needed to plot efficiency as a function of the given bin variable from list of EffInfos'''
  graphValues: List[Tuple[float, UFloat]] = [(effInfo.binInfo.centers[binVarName], effInfo.value)
    for effInfo in effInfos if binVarName in effInfo.binInfo.varNames]
  return graphValues


def plotEfficiencies1D(
  efficiencies:      List[EffInfo],
  binningVar:        str,  # name of binning variable to plot
  pdfDirName:        str,  # directory name the PDF file will be written to
  pdfFileNameSuffix: str   = "",
  particle:          str   = "Proton",
  channel:           str   = "4pi",
  markerSize:        float = 0.75,
) -> None:
  '''Plots efficiency as a function of `binningVar` for 1-dimensional binning'''
  print(f"Plotting efficiency as a function of binning variable '{binningVar}'")
  efficiencyGraph = plotFitResults.getParValueGraph1D(getEffValuesForGraph1D(binningVar, efficiencies))
  efficiencyGraph.SetTitle(f"{particle} Track-Finding Efficiency ({channel})")
  efficiencyGraph.SetMarkerStyle(ROOT.kFullCircle)  # type: ignore
  efficiencyGraph.SetMarkerSize(markerSize)
  assert binningVar in BINNING_VAR_PLOT_INFO, f"No plot information for binning variable '{binningVar}'"
  efficiencyGraph.GetXaxis().SetTitle(f"{BINNING_VAR_PLOT_INFO[binningVar]['label']} ({BINNING_VAR_PLOT_INFO[binningVar]['unit']})")
  efficiencyGraph.GetYaxis().SetTitle("Efficiency")
  efficiencyGraph.SetMinimum(0)
  efficiencyGraph.SetMaximum(1)
  canv = ROOT.TCanvas(f"{particle}_{channel}_mm2_eff_{binningVar}{pdfFileNameSuffix}", "")  # type: ignore
  efficiencyGraph.Draw("AP")
  # #TODO? # indicate value from fit of overall distributions
  # line = ROOT.TLine()  # type: ignore
  # line.SetLineStyle(ROOT.kDashed)  # type: ignore
  # line.DrawLine(efficienciesKinBinsGraph.GetXaxis().GetXmin(), overallEff.nominal_value, efficienciesKinBinsGraph.GetXaxis().GetXmax(), overallEff.nominal_value)
  # # indicate weighted average of efficiencies in kinematic bins
  # meanEff = np.average(yVals, weights = [1 / (yErr**2) for yErr in yErrs])
  # line.SetLineColor(ROOT.kRed + 1)  # type: ignore
  # line.DrawLine(efficiencyGraph.GetXaxis().GetXmin(), meanEff, efficiencyGraph.GetXaxis().GetXmax(), meanEff)
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)  # type: ignore
  makePlots.setupPlotStyle()
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  # type: ignore  #TODO use BRUFIT environment variable

  # echo and parse command line
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Plots BruFit results.")
  parser.add_argument("outputDirName", type = str, nargs = "?", default = "./BruFitOutput", help = "The path to the BruFit output directory; (default: '%(default)s')")
  args = parser.parse_args()
  dataSets    = ["Total", "Found", "Missing"]
  fitVariable = "MissingMassSquared_Measured"  #TODO get this info from BruFit's ROOT.Setup

  for readIntegrals in (True, False):
    print("Calculating efficiencies from " + ("integrals of data histograms" if readIntegrals else "fit results"))
    yieldInfos:  Dict[str, List[ParInfo]]           = {}    # parInfos[<dataset>][<bin>]
    binVarNames: Union[List[Tuple[str, ...]], None] = None  # binning variables for each binning
    for dataSet in dataSets:
      print("Reading " + ("integrals of data histograms" if readIntegrals else "yields") + f" for '{dataSet}' dataset")
      fitResultDirName = f"{args.outputDirName}/{dataSet}"
      yieldInfos[dataSet] = readYieldInfosForBinning(BinningInfo([], fitResultDirName), fitVariable, readIntegrals)  # overall yield values
      binVarNamesInDataSet = []
      for binningInfo in plotFitResults.getBinningInfosFromDir(fitResultDirName):
        if binningInfo:
          binVarNamesInDataSet.append(binningInfo.varNames)
          yieldInfos[dataSet][len(yieldInfos[dataSet]):] = readYieldInfosForBinning(binningInfo, fitVariable, readIntegrals)  # append yield values
      if binVarNames is None:
        binVarNames = binVarNamesInDataSet
      else:
        assert binVarNamesInDataSet == binVarNames, f"The binning variables {binVarNamesInDataSet} of dataset '{dataSet}' are different from the binning variables {binVarNames} of the previous one"
    effInfos = calculateEfficiencies(yieldInfos)
    if effInfos:
      if effInfos[0].binInfo.name == "Overall":
        print("Overall efficiency from " + ("data-histogram integrals" if readIntegrals else "fits") + f" is {effInfos[0].value}")
      if binVarNames:
        for binningVars in binVarNames:
          if len(binningVars) == 1:
            plotEfficiencies1D(effInfos, binningVars[0], args.outputDirName, "_integral" if readIntegrals else "")
