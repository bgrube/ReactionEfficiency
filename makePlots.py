#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetHistMinimumZero(True)
ROOT.gROOT.LoadMacro("~/rootlogon.C")


def getCanvName(fileName, histName):
  histName = histName.replace("MissingMassSquared/", "")
  particle = "Pi+" if "pipmiss" in fileName else ('Pi-' if "pimmiss" in fileName else "Proton")
  channel  = "4pi" if fileName.count("pi") == 4 else "2pi"
  dataset  = fileName.split('.')[1]
  if "missingmass" in histName.lower():
    return f"justin_{particle}_{channel}_{'mm2' if 'Squared' in histName else 'mm'}{'_sb' if histName.endswith('Sideband') else ''}_{dataset}"
  else:
    return f"justin_{particle}_{channel}_{histName.lower()}_{dataset}"


# missing particle kinematics
histNames = [
  "MissingParticleMomVsTheta",
  "MissingParticlePhiVsTheta"]

inFiles = [ROOT.TFile(inFileName) for inFileName in inFileNames]
hists = [[inFile.Get(histName) for histName in histNames] for inFile in inFiles]
canvNames = [[getCanvName(inFileName, histName) for histName in histNames] for inFileName in inFileNames]
print(canvNames)
canvs = [[ROOT.TCanvas(name, name, 600, 600) for name in names] for names in canvNames]
print(canvs)
for fileIndex, histsInFile in enumerate(hists):
  for histIndex, hist in enumerate(histsInFile):
    canv = canvs[fileIndex][histIndex]
    canv.cd()
    hist.SetStats(False)
    hist.Draw("COLZ")
    hist.GetXaxis().SetRangeUser(0, 90)
    canv.SaveAs(".pdf")
