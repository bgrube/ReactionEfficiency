#!/usr/bin/env python3


from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
import copy
from dataclasses import replace
import functools
import numpy as np
import os
from scipy import stats

from uncertainties import UFloat, ufloat

import ROOT
if __name__ == "__main__":
  ROOT.PyConfig.DisableRootLogon = True  # do not change style of canvases loaded from fit result files

from plotFitResults import (
  acceptFitResult,
  BinInfo,
  BinningInfo,
  getAxisInfoForBinningVar,
  getBinningInfosFromDir,
  ParInfo,
)
from plotTools import (
  printGitInfo,
  redrawFrame,
  setupPlotStyle,
)


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


#TODO move to and use in plotFitResults.py
def readFitParameters(
  fitResultFileName: str,
  binInfo:           BinInfo,
) -> tuple[ParInfo | None, bool]:
  """Reads fit parameter values from a fit-result file"""
  if not os.path.isfile(fitResultFileName):
    print(f"Cannot find file '{fitResultFileName}'. Skipping.")
    return None, False
  fitResultObjName = "MinuitResult"
  print(f"Reading fit result object '{fitResultObjName}' from file '{fitResultFileName}'")
  fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
  fitResult = fitResultFile.Get(fitResultObjName)
  # read all parameters
  parValuesInBin: dict[str, UFloat] = {fitPar.GetName(): ufloat(fitPar.getVal(), fitPar.getError()) for fitPar in fitResult.floatParsFinal()}
  fitResultFile.Close()
  return ParInfo(binInfo, fitVariable, parValuesInBin), acceptFitResult(fitResult)


def plotBootstrapDistribution(
  parName:        str,                # name of fit parameter
  parInfosBs:     Sequence[ParInfo],  # bootstrap samples of fit parameters
  parInfoNominal: ParInfo | None,     # nominal parameter estimates
  dataSet:        str,                # name of dataset
  nmbBins:        int = 20,           # number of histogram bins
) -> ParInfo | None:  # ParInfo with deviation of bootstrap estimate from nominal estimate
  """Plots bootstrap distribution of a fit parameter"""
  print(f"Plotting bootstrap distribution for parameter '{parName}' and dataset '{dataSet}'")
  parValuesBs = np.array([parInfoBs.values[parName].nominal_value for parInfoBs in parInfosBs], dtype = np.float64)
  min = np.min(parValuesBs)
  max = np.max(parValuesBs)
  halfRange = (max - min) * 1.1 / 2.0
  center = (min + max) / 2.0
  centers: dict[str, float] = dict(parInfosBs[0].binInfo.centers)
  histBs = ROOT.TH1D(
    f"bootstrap_{dataSet}_"
    + (("_".join((str(item) for center in centers.items() for item in center)) + "_") if centers else "")
    + f"{parName}",
    f"{dataSet}" + (f", {centers}" if centers else "") + f";{parName};Count",
    nmbBins, center - halfRange, center + halfRange
  )
  # fill histogram
  np.vectorize(histBs.Fill, otypes = [int])(parValuesBs)
  # draw histogram
  canv = ROOT.TCanvas()
  histBs.SetMinimum(0)
  histBs.SetLineColor(ROOT.kBlue + 1)
  histBs.Draw("E")
  # indicate bootstrap estimate
  meanBs   = np.mean(parValuesBs)
  stdDevBs = np.std(parValuesBs, ddof = 1)
  yCoord = histBs.GetMaximum() / 4
  markerBs = ROOT.TMarker(meanBs, yCoord, ROOT.kFullCircle)
  markerBs.SetMarkerColor(ROOT.kBlue + 1)
  markerBs.SetMarkerSize(0.75)
  markerBs.Draw()
  lineBs = ROOT.TLine(meanBs - stdDevBs, yCoord, meanBs + stdDevBs, yCoord)
  lineBs.SetLineColor(ROOT.kBlue + 1)
  lineBs.Draw()
  meanEst = stdDevEst = None
  if parInfoNominal is not None and parInfoNominal.values:
    # indicate nominal estimate
    meanEst   = parInfoNominal.values[parName].nominal_value
    stdDevEst = parInfoNominal.values[parName].std_dev
    markerEst = ROOT.TMarker(meanEst,  yCoord / 2, ROOT.kFullCircle)
    markerEst.SetMarkerColor(ROOT.kGreen + 2)
    markerEst.SetMarkerSize(0.75)
    markerEst.Draw()
    lineEst = ROOT.TLine(meanEst - stdDevEst, yCoord / 2, meanEst + stdDevEst, yCoord / 2)
    lineEst.SetLineColor(ROOT.kGreen + 2)
    lineEst.Draw()
    # plot Gaussian that corresponds to estimate from uncertainty propagation
    gaussian = ROOT.TF1("gaussian", "gausn(0)", center - halfRange, center + halfRange)
    gaussian.SetParameters(len(parValuesBs) * histBs.GetBinWidth(1), meanEst, stdDevEst)
    gaussian.SetLineColor(ROOT.kGreen + 2)
    gaussian.Draw("SAME")
    # print chi^2 of Gaussian and histogram
    chi2     = histBs.Chisquare(gaussian, "L")
    chi2Prob = stats.distributions.chi2.sf(chi2, nmbBins)
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    label.SetTextSize(0.04)
    label.SetTextColor(ROOT.kGreen + 2)
    label.DrawLatex(0.13, 0.85, f"#it{{#chi}}^{{2}}/n.d.f. = {chi2:.2f}/{nmbBins}, prob = {chi2Prob * 100:.0f}%")
  # add legend
  legend = ROOT.TLegend(0.7, 0.75, 0.99, 0.95)
  legend.AddEntry(histBs, "Bootstrap samples", "LE")
  entry = legend.AddEntry(markerBs, "Bootstrap estimate", "LP")
  entry.SetLineColor(ROOT.kBlue + 1)
  legend.AddEntry(0, f"Mean = {meanBs}",      "")
  legend.AddEntry(0, f"Std. dev. = {stdDevBs}", "")
  if parInfoNominal is not None and parInfoNominal.values:
    entry = legend.AddEntry(markerEst, "Nominal estimate", "LP")
    entry.SetLineColor(ROOT.kGreen + 2)
    legend.AddEntry(gaussian, "Nominal estimate Gaussian", "LP")
  legend.Draw()
  canv.SaveAs(f"{parInfosBs[0].binInfo.dirName}/{histBs.GetName()}.pdf")
  if parInfoNominal is None or meanEst is None or stdDevEst is None:
    return None
  # return relative difference of nominal estimate from bootstrap estimate
  parInfoDeviation = replace(
    parInfoNominal,
    values = {
      f"{parName}_meanDev"        : ufloat((meanBs   - meanEst  ) / stdDevBs, 0),
      f"{parName}_uncertaintyDev" : ufloat((stdDevBs - stdDevEst) / stdDevBs, 0),
    },
  )
  return parInfoDeviation


if __name__ == "__main__":
  printGitInfo()
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")
  setupPlotStyle()
  ROOT.gStyle.SetOptStat(False)
  ROOT.gStyle.SetPalette(ROOT.kLightTemperature)


  nmbBootstrapSamples    = 100  #TODO determine from files
  outputDirNameBs        = "./fits.Bs_100/2017_01-ver03/BruFitOutput.data_2017_01-ver03_allFixed"  # fit directory with bootstrap samples
  outputDirNameNominal   = "./fits.nominal/2017_01-ver03/BruFitOutput.data_2017_01-ver03_allFixed"  # fit directory with nominal fit results
  dataSets               = ["Total", "Found", "Missing"]
  fitVariable            = "MissingMassSquared_Measured"
  acceptedFitResultsOnly = False

  for dataSet in dataSets:  # loop over datasets
    fitResultDirName  = f"{outputDirNameBs}/{dataSet}"
    binningInfoOverall = BinningInfo(
      infos   = [BinInfo(name = "", binDefs = {}, dirName = fitResultDirName)],
      dirName = fitResultDirName,
    )  # dummy binning info for overall distribution
    for binningInfo in [binningInfoOverall] + getBinningInfosFromDir(fitResultDirName):  # loop over kinematic binnings
      if binningInfo:
        print(f"Reading bootstrap distribution for '{dataSet}' dataset and binning {binningInfo}")
        deviations: defaultdict[str, list[ParInfo]] = defaultdict(list)  # parameter deviations for all bins [parName][binIndex]
        for binInfoBs in binningInfo:  # loop over kinematic bins
          # read nominal parameter estimates
          binInfoNominal = copy.deepcopy(binInfoBs)
          binInfoNominal.dirName = binInfoBs.dirName.replace(outputDirNameBs, outputDirNameNominal, 1)
          print(f"Reading nominal estimate for bin {binInfoNominal}")
          parInfoNominal, acceptedFitResult = readFitParameters(binInfoNominal.fitResultFileName, binInfoNominal)
          if parInfoNominal is None or not parInfoNominal.values:
            print(f"No nominal parameter estimates found for dataset '{dataSet}' and bin {binInfoBs}. Skipping.")
            continue
          # read bootstrap samples
          binInfoBs.nmbBootstrapSamples = nmbBootstrapSamples
          print(f"Reading bootstrap distribution for bin {binInfoBs}")
          parInfosBs: list[ParInfo] = []  # parameter values for all bootstrap samples
          for bootstrapIndex, fitResultFileName in binInfoBs.bootstrapFileNames:  # loop over bootstrap samples
            parInfoBs, _  = readFitParameters(fitResultFileName, binInfoBs)
            if parInfoBs is not None and parInfoBs.values:
              print(f"Parameter values for '{dataSet}' dataset, bin {binInfoBs}, and bootstrap index {bootstrapIndex}: {parInfoBs}")
              parInfosBs.append(parInfoBs)
          if not parInfosBs:
            print(f"No bootstrap samples found for dataset '{dataSet}' and bin {binInfoBs}. Skipping.")
            continue
          # plot bootstrap distribution
          parNames = tuple(parInfoNominal.values.keys())  # assume parInfoNominal and parInfosBs have the same parameters
          for parName in parNames:
            parInfoDeviation = plotBootstrapDistribution(parName, parInfosBs, parInfoNominal, dataSet)
            if parInfoDeviation is not None and not (acceptedFitResultsOnly and not acceptedFitResult):
              deviations[parName].append(parInfoDeviation)
        # plot deviations of nominal estimates from bootstrap estimates for all parameters
        if len(binningInfo.varNames) == 2:
          binningVars = binningInfo.varNames[:2]
          for parName, deviationsForParName in deviations.items():  # loop over parameters
            for devName in ("mean", "uncertainty"):
              canv = ROOT.TCanvas()
              binningVarLabels: list[str] = [""] * len(binningVars)
              binningVarUnits:  list[str] = [""] * len(binningVars)
              for index, binningVar in enumerate(binningVars):
                _, binningVarLabels[index], binningVarUnits[index] = getAxisInfoForBinningVar(binningVar)
              hist = ROOT.TH2D(
                f"bootstrap_{dataSet}_{parName}_{devName}Dev",
                f"{dataSet} {parName} {devName} Relative Difference"
                f";{binningVarLabels[0]} ({binningVarUnits[0]})"
                f";{binningVarLabels[1]} ({binningVarUnits[1]})"
                f";({devName}_{{BS}} #minus {devName}_{{Fit}}) / #sigma_{{BS}}",
                binningInfo.varNmbBins(binningVars[0]), *binningInfo.varRange(binningVars[0]),
                binningInfo.varNmbBins(binningVars[1]), *binningInfo.varRange(binningVars[1]),
              )
              # fill histogram
              for deviation in deviationsForParName:
                hist.SetBinContent(
                  hist.FindBin(deviation.binInfo.center(binningVars[0]), deviation.binInfo.center(binningVars[1])),
                  deviation.values[f"{parName}_{devName}Dev"].nominal_value,
                )
              # draw histogram
              hist.SetMinimum(-0.5)
              hist.SetMaximum(+0.5)
              ROOT.gStyle.SetPaintTextFormat("1.3f")
              hist.Draw("COLZ TEXT")
              hist.SetStats(False)
              canv.SetRightMargin(0.15)
              redrawFrame(canv)
              canv.SaveAs(f"{binningInfo.dirName}/{hist.GetName()}.pdf")
