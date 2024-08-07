#!/usr/bin/env python3


from __future__ import annotations

from collections import defaultdict
import functools
import numpy as np
import nptyping as npt
import os

from uncertainties import UFloat, ufloat

import ROOT
if __name__ == "__main__":
  ROOT.PyConfig.DisableRootLogon = True  # do not change style of canvases loaded from fit result files

from plotFitResults import (
  BinInfo,
  ParInfo,
)
from plotTools import (
  printGitInfo,
  setupPlotStyle,
)


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


def plotBootstrapDistribution(
  parName:       str,  # name of fit parameter
  parValues:     npt.NDArray[npt.Shape["*, *"], npt.Float64],  # array with bootstrap samples of fit parameter
  dataSet:       str,  # name of dataset
  outputDirName: str,  # directory where pdf is written to
  nmbBins:       int = 20,  # number of histogram bins
):
  """Plots bootstrap distribution of a fit parameter"""
  print(f"Plotting bootstrap distribution for parameter '{parName}' and dataset '{dataSet}'")
  min = np.min(parValues)
  max = np.max(parValues)
  halfRange = (max - min) * 1.1 / 2.0
  center = (min + max) / 2.0
  histBs = ROOT.TH1D(f"bootstrap_{dataSet}_{parName}", f"{dataSet};{parName};Count",
                     nmbBins, center - halfRange, center + halfRange)
  # fill histogram
  np.vectorize(histBs.Fill, otypes = [int])(parValues)
  # draw histogram
  canv = ROOT.TCanvas()
  histBs.SetMinimum(0)
  histBs.SetLineColor(ROOT.kBlue + 1)
  histBs.Draw("E")
  # indicate bootstrap estimate
  meanBs   = np.mean(parValues)
  stdDevBs = np.std(parValues, ddof = 1)
  yCoord = histBs.GetMaximum() / 4
  markerBs = ROOT.TMarker(meanBs, yCoord, ROOT.kFullCircle)
  markerBs.SetMarkerColor(ROOT.kBlue + 1)
  markerBs.SetMarkerSize(0.75)
  markerBs.Draw()
  lineBs = ROOT.TLine(meanBs - stdDevBs, yCoord, meanBs + stdDevBs, yCoord)
  lineBs.SetLineColor(ROOT.kBlue + 1)
  lineBs.Draw()
  # add legend
  legend = ROOT.TLegend(0.7, 0.8, 0.99, 0.99)
  legend.AddEntry(histBs, "Bootstrap samples", "LE")
  entry = legend.AddEntry(markerBs, "Bootstrap estimate", "LP")
  entry.SetLineColor(ROOT.kBlue + 1)
  legend.AddEntry(0, f"Mean = {meanBs}",      "")
  legend.AddEntry(0, f"Uncert. = {stdDevBs}", "")
  legend.Draw()
  canv.SaveAs(f"{outputDirName}/{dataSet}/{histBs.GetName()}.pdf")


if __name__ == "__main__":
  printGitInfo()
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")
  setupPlotStyle()
  ROOT.gStyle.SetOptStat(False)

  outputDirName = "/home/bgrube/Analysis/ProtonTrackEfficiency/ReactionEfficiency/fits.new/2017_01-ver03_goodToF/noShowers/BruFitOutput.data_2017_01-ver03_goodToF_allFixed"
  nmbBootstrapSamples = 10
  dataSets = ["Total", "Found", "Missing"]
  fitVariable = "MissingMassSquared_Measured"

    # plot fits and read parameter values from fit results
  for dataSet in dataSets:
    fitResultDirName  = f"{outputDirName}/{dataSet}"
    print(f"Reading bootstrap distribution for overall fit result for '{dataSet}' dataset")
    binInfo = BinInfo("", {}, {}, fitResultDirName, nmbBootstrapSamples)
    parInfos: list[ParInfo] = []  # parameter values for all bootstrap samples
    for bootstrapIndex, fitResultFileName in binInfo.bootstrapFileNames:
      if not os.path.isfile(fitResultFileName):
        print(f"Cannot find file '{fitResultFileName}'. Skipping bin {binInfo}.")
        continue
      print(f"Reading fit result object 'MinuitResult' from file '{fitResultFileName}'")
      fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
      fitResult     = fitResultFile.Get("MinuitResult")
      # read all parameters
      parValuesInBin: dict[str, UFloat] = {fitPar.GetName() : ufloat(fitPar.getVal(), fitPar.getError()) for fitPar in fitResult.floatParsFinal()}
      fitResultFile.Close()
      parInfo = ParInfo(binInfo, fitVariable, parValuesInBin)
      if parInfo is not None and parInfo.values:
        print(f"Parameter values for bootstrap index {bootstrapIndex}: {parInfo}")
        parInfos.append(parInfo)
      if not parInfos:
        print(f"No parameter values found for dataset '{dataSet}'")
        continue
    # plot bootstrap distributions
    parNames = tuple(parInfos[0].values.keys())
    for parName in parNames:
      parValues = np.array([parInfo.values[parName].nominal_value for parInfo in parInfos], dtype = np.float64)
      plotBootstrapDistribution(
        parName,
        parValues,
        dataSet,
        outputDirName,
      )
