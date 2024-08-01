#!/usr/bin/env python3

import numpy as np
import ROOT

ROOT.PyConfig.DisableRootLogon = True

if __name__ == "__main__":
  # create TH3 to serve as a frame
  canv = ROOT.TCanvas()
  frame = ROOT.TH3F("frame", ";x;y;z", 10, 0, 10, 10, 0, 10, 10, 0, 10)
  frame.SetStats(False)
  frame.Draw()

  # create TGraph2D and draw it into TH3 frame
  xVals = np.array([1, 3, 5, 7, 9], dtype = "d")
  yVals = np.array([1, 9, 1, 9, 1], dtype = "d")
  zVals = np.array([3, 7, 3, 7, 3], dtype = "d")
  graph = ROOT.TGraph2D(len(xVals), xVals, yVals, zVals)
  graph.SetMarkerStyle(20)
  graph.Draw("P0 SAME")
  # canv.Update()  # <- if this line is commented out, the bug does not occur
  canv.SaveAs("TH3Graph2DBug.png")
