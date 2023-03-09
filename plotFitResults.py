#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import os

import ROOT
ROOT.PyConfig.DisableRootLogon = True

import plotEfficiencies


def plotFitResults(fitResultDirName):
  bins = plotEfficiencies.getBinningFromFile(fitResultDirName)
  binNames = [str(binName) for binName in bins.GetBinNames()]
  fitResultFileNames = [str(fileName).replace("TreeData.root", "ResultsHSMinuit2.root") for fileName in bins.GetFileNames()]
  for index, fitResultFileName  in enumerate(fitResultFileNames):
    fitResultFile = ROOT.TFile.Open(fitResultFileName, "READ")
    print(f"Plotting fit result in '{fitResultFileName}'")
    canvName = f"{binNames[index]}_{fitVariable}"
    canv = fitResultFile.Get(canvName)
    # remove TPaveText with fit parameters
    dataFitPad = canv.GetListOfPrimitives().FindObject(f"{canvName}_1")
    paramBox = dataFitPad.GetListOfPrimitives().FindObject(f"{binNames[index]}TotalPDF_paramBox")
    dataFitPad.GetListOfPrimitives().Remove(paramBox)
    # paramBox.SetBorderSize(0)
    # paramBox.SetFillStyle(0)
    canv.SaveAs(f"{fitResultDirName}/{canv.GetName()}.pdf")
    fitResultFile.Close()


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable

  outputDirName = "BruFitOutput"
  dataSets      = ["Total", "Found", "Missing"]
  fitVariable   = "MissingMassSquared_Measured"  #TODO read this from ROOT.Setup

  for dataSet in dataSets:
    # plot fit results for kinematic bins
    plotFitResults(f"{outputDirName}/{dataSet}")
