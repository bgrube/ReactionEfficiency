#!/usr/bin/env python3


from __future__ import annotations

from collections.abc import Sequence
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
  BinningInfo,
  getBinningInfosFromDir,
  ParInfo,
)
from plotTools import (
  printGitInfo,
  setupPlotStyle,
)


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


def plotBootstrapDistribution(
  parName:  str,  # name of fit parameter
  parInfos: Sequence[ParInfo],  # bootstrap samples of fit parameters
  dataSet:  str,  # name of dataset
  nmbBins:  int = 20,  # number of histogram bins
):
  """Plots bootstrap distribution of a fit parameter"""
  print(f"Plotting bootstrap distribution for parameter '{parName}' and dataset '{dataSet}'")
  parValues = np.array([parInfo.values[parName].nominal_value for parInfo in parInfos], dtype = np.float64)
  min = np.min(parValues)
  max = np.max(parValues)
  halfRange = (max - min) * 1.1 / 2.0
  center = (min + max) / 2.0
  centers: dict[str, float] = parInfos[0].binInfo.centers
  histBs = ROOT.TH1D(f"bootstrap_{dataSet}_"
                     + (("_".join((str(item) for center in centers.items() for item in center)) + "_") if centers else "")
                     + f"{parName}",
                     f"{dataSet}" + (f", {centers}" if centers else "") + f";{parName};Count",
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
  canv.SaveAs(f"{parInfos[0].binInfo.dirName}/{histBs.GetName()}.pdf")


if __name__ == "__main__":
  printGitInfo()
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")
  setupPlotStyle()
  ROOT.gStyle.SetOptStat(False)

  outputDirName = "/home/bgrube/Analysis/ProtonTrackEfficiency/ReactionEfficiency/fits.Bs_10/2017_01-ver03_goodToF/noShowers/BruFitOutput.data_2017_01-ver03_goodToF_allFixed"
  nmbBootstrapSamples = 10  #TODO determine from files
  dataSets = ["Total", "Found", "Missing"]
  fitVariable = "MissingMassSquared_Measured"

  for dataSet in dataSets:  # loop over datasets
    fitResultDirName  = f"{outputDirName}/{dataSet}"
    binningInfoOverall = BinningInfo(
      infos   = [BinInfo(name = "", centers = {}, widths = {}, dirName = fitResultDirName)],
      dirName = fitResultDirName,
    )  # dummy binning info for overall distribution
    for binningInfo in [binningInfoOverall] + getBinningInfosFromDir(fitResultDirName):  # loop over kinematic binnings
      if binningInfo:
        print(f"Reading bootstrap distribution for '{dataSet}' dataset and binning {binningInfo}")
        for binInfo in binningInfo:  # loop over kinematic bins
          binInfo.nmbBootstrapSamples = nmbBootstrapSamples
          print(f"Reading bootstrap distribution for bin {binInfo}")
          parInfos: list[ParInfo] = []  # parameter values for all bootstrap samples
          for bootstrapIndex, fitResultFileName in binInfo.bootstrapFileNames:  # loop over bootstrap samples
            if not os.path.isfile(fitResultFileName):
              print(f"Cannot find file '{fitResultFileName}'. Skipping bootstrap sample {bootstrapIndex}.")
              continue
            print(f"Reading fit result object 'MinuitResult' from file '{fitResultFileName}'")
            fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
            fitResult     = fitResultFile.Get("MinuitResult")
            # read all parameters
            parValuesInBin: dict[str, UFloat] = {fitPar.GetName() : ufloat(fitPar.getVal(), fitPar.getError()) for fitPar in fitResult.floatParsFinal()}
            fitResultFile.Close()
            parInfo = ParInfo(binInfo, fitVariable, parValuesInBin)
            if parInfo is not None and parInfo.values:
              print(f"Parameter values for '{dataSet}' dataset, bin {binInfo}, and bootstrap index {bootstrapIndex}: {parInfo}")
              parInfos.append(parInfo)
          if not parInfos:
            print(f"No parameter values found for dataset '{dataSet}' and bin {binInfo}")
            continue
          # plot bootstrap distributions
          parNames = tuple(parInfos[0].values.keys())
          for parName in parNames:
            plotBootstrapDistribution(parName, parInfos, dataSet)
