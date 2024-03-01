#!/usr/bin/env python3


import functools

import ROOT

import plotTools


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


if __name__ == "__main__":
  fitDirName  = "./out"
  fitVariable = "MissingMassSquared_Measured"  # name of branch that holds data to fit and template-data for signal and background, respectively

  print(">>> plotting fit result")
  fitResultFile = ROOT.TFile.Open(f"{fitDirName}/ResultsHSMinuit2.root", "READ")
  canvName = f"_{fitVariable}"
  canv = fitResultFile.Get(canvName)
  # improve TPaveText with fit parameters
  dataFitPad = canv.GetListOfPrimitives().FindObject(f"{canvName}_1")
  paramBox = dataFitPad.GetListOfPrimitives().FindObject("TotalPDF_paramBox")
  # remove filled frame
  paramBox.SetBorderSize(0)
  paramBox.SetFillStyle(0)
  plotTools.drawZeroLine(dataFitPad)
  canv.SaveAs(f"{fitDirName}/ResultsHSMinuit2.pdf")
  fitResultFile.Close()
