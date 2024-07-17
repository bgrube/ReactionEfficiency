#!/usr/bin/env python3

import ROOT
import subprocess

ROOT.PyConfig.DisableRootLogon = True
ROOT.gROOT.SetBatch()

def createTestFile():
  myColor = ROOT.TColor(ROOT.TColor.GetFreeColorIndex(), 1, 1, 1, "myColor", 1)
  myColor.Print()
  myColorIndex = myColor.GetNumber()
  ROOT.gStyle.SetFillColor(myColorIndex)
  ROOT.gStyle.SetPalette(ROOT.kViridis)  # <- if this line is commented out, the bug does not occur
  canv = ROOT.TCanvas("canv", "canv")
  canv.SaveAs("TColorBug.root")

def readTestFile():
  file = ROOT.TFile.Open("TColorBug.root", "READ")
  file.Get("canv")
  freeColorIndex = ROOT.TColor.GetFreeColorIndex()
  print(f"ROOT.TColor.GetFreeColorIndex() = {freeColorIndex}")
  try:
    ROOT.gROOT.GetColor(freeColorIndex).Print()
  except:
    print(f"Expected behavior: index {freeColorIndex} returned by ROOT.TColor.GetFreeColorIndex() is not used")
  else:
    print(f"Wrong behavior: index {freeColorIndex} returned by ROOT.TColor.GetFreeColorIndex() is already used")

if __name__ == "__main__":
  print("Creating test file in a separate process")
  subprocess.run("python3 -c 'import TColorBug; TColorBug.createTestFile()'", shell = True)
  print("Reading test file")
  readTestFile()
