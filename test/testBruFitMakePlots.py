#!/usr/bin/env python3


import argparse
import functools

import ROOT

import plotTools


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Plots BruFit test result")
  parser.add_argument("fitDirName",  type = str, help = "The path to the BruFit output directory")
  parser.add_argument("fitVariable", type = str, help = "The name of branch that holds data to fit")
  args = parser.parse_args()

  fitResultFile = ROOT.TFile.Open(f"{args.fitDirName}/ResultsHSMinuit2.root", "READ")
  canvName = f"_{args.fitVariable}"
  canv = fitResultFile.Get(canvName)
  # improve TPaveText with fit parameters
  dataFitPad = canv.GetListOfPrimitives().FindObject(f"{canvName}_1")
  paramBox = dataFitPad.GetListOfPrimitives().FindObject("TotalPDF_paramBox")
  # remove filled frame
  paramBox.SetBorderSize(0)
  paramBox.SetFillStyle(0)
  plotTools.drawZeroLine(dataFitPad)
  canv.SaveAs(f"{args.fitDirName}/ResultsHSMinuit2.pdf")
  fitResultFile.Close()
