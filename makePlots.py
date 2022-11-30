#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetHistMinimumZero(True)
ROOT.gROOT.LoadMacro("~/rootlogon.C")


def getCanvName(fileName, histName):
  particle = "Pi+" if "pipmiss" in fileName else ('Pi-' if "pimmiss" in fileName else "Proton")
  channel  = "4pi" if fileName.count("pi") == 4 else "2pi"
  dataset  = fileName.split('.')[1]
  if "missingmass" in histName.lower():
    return f"justin_{particle}_{channel}_{'mm2' if 'Squared' in histName else 'mm'}{'_sb' if histName.endswith('Sideband') else ''}_{dataset}"
  else:
    return f"justin_{particle}_{channel}_{histName.lower()}_{dataset}"


# total missing-mass (squared) distributions
inFileNames = [
  "pippimpmiss.30370.root",
  "pippimpmiss.30370_acc.root"]
  # "pippippimpimpmiss.30370.root",
  # "pippippimpimpmiss.30370_acc.root"]
histNames = [
  "MissingMass",
  "MissingMassSideband",
  "MissingMassSquared",
  "MissingMassSquaredSideband"]
rebinFactor = 20

inFiles = [ROOT.TFile(inFileName) for inFileName in inFileNames]
hists = [[inFile.Get(histName) for histName in histNames] for inFile in inFiles]
canvNames = [[getCanvName(inFileName, histName) for histName in histNames] for inFileName in inFileNames]
canvs = [[ROOT.TCanvas(name, name, 1470, 891) for name in names] for names in canvNames]
for fileIndex, histsInFile in enumerate(hists):
  for histIndex, hist in enumerate(histsInFile):
    canv = canvs[fileIndex][histIndex]
    canv.cd()
    if "Squared" in hist.GetName():
      mmRange = (-0.5, 4)
    else:
      mmRange = (0, 2)
    print(f"{canv.GetName()}: {hist.Integral(hist.FindBin(mmRange[0]), hist.FindBin(mmRange[1]))} events in range {str(mmRange)}")
    hist.Rebin(rebinFactor)
    hist.GetXaxis().SetRangeUser(*mmRange)
    hist.Draw("HIST")
    if histIndex % 2 == 0:
      # overlay RF sidebands
      sbHist = hists[fileIndex][histIndex + 1]
      if sbHist.Integral() != 0:
        sbHistCopy = sbHist.DrawCopy("HIST SAME")
        sbHistCopy.Rebin(rebinFactor)
        sbHistCopy.SetLineColor(ROOT.kGreen+2)
        if sbHistCopy.GetMaximum() > hist.GetMaximum():
          hist.SetMaximum(1.1 * sbHistCopy.GetMaximum())
    canv.SaveAs(".pdf")
    print()

inFileNames = [
  "pippimpmiss.30370_acc.root",
  "pippippimpimpmiss.30370_acc.root"]
histName = "MissingMassSquared"

inFiles = [ROOT.TFile(inFileName) for inFileName in inFileNames]
hists = [inFile.Get(histName) for inFile in inFiles]
canvNames = [f"{getCanvName(inFileName, histName)}_2" for inFileName in inFileNames]
canvs = [ROOT.TCanvas(canvName, canvName, 1470, 891) for canvName in canvNames]
for index, hist in enumerate(hists):
  canv = canvs[index]
  canv.cd()
  hist.Rebin(rebinFactor)
  hist.GetXaxis().SetRangeUser(-0.5, 4)
  hist.Draw("HIST")
  canv.SaveAs(".pdf")

# RF weights
histName = "RFWeight"

inFiles = [ROOT.TFile(inFileName) for inFileName in inFileNames]
hists = [inFile.Get(histName) for inFile in inFiles]
canvNames = [f"{getCanvName(inFileName, histName)}" for inFileName in inFileNames]
canvs = [ROOT.TCanvas(canvName, canvName, 1470, 891) for canvName in canvNames]
for index, hist in enumerate(hists):
  canv = canvs[index]
  canv.cd()
  canv.SetLogy()
  hist.Draw("HIST")
  canv.Update()  # needed otherwise TPaveStats object is not created
  stats = hist.FindObject("stats")
  # force class
  # stats.__class__ = ROOT.TPaveStats
  stats.SetOptStat(111111)
  canv.SaveAs(".pdf")


# missing particle kinematics
histNames = [
  "MissingParticleMomVsTheta",
  "MissingParticlePhiVsTheta"]

inFiles = [ROOT.TFile(inFileName) for inFileName in inFileNames]
hists = [[inFile.Get(histName) for histName in histNames] for inFile in inFiles]
canvNames = [[getCanvName(inFileName, histName) for histName in histNames] for inFileName in inFileNames]
print(canvNames)
canvs = [[ROOT.TCanvas(name, name, 1470, 891) for name in names] for names in canvNames]
print(canvs)
for fileIndex, histsInFile in enumerate(hists):
  for histIndex, hist in enumerate(histsInFile):
    canv = canvs[fileIndex][histIndex]
    canv.cd()
    hist.SetStats(False)
    hist.Draw("COLZ")
    hist.GetXaxis().SetRangeUser(0, 90)
    canv.SaveAs(".pdf")
