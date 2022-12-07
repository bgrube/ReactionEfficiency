#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetHistMinimumZero(True)
ROOT.gROOT.LoadMacro("~/rootlogon.C")


def plotTopologies(maxNmbTopologies = 10, normalize = False):
  inFileName = "pippippimpimpmiss.root"
  histBaseName = "MissingMassSquared/ThrownTopologies"
  colors = {
    "Found"   : ROOT.kGreen + 2,
    "Missing" : ROOT.kRed + 1,
    ""        : ROOT.kBlue
  }

  inFile = ROOT.TFile(inFileName)
  cases = ["Found", "Missing", ""]
  hists = {case : inFile.Get(histBaseName + ("_" + case if case != "" else "")) for case in cases}
  if normalize:
    for hist in hists.values():
      hist.Scale(100 / hist.Integral())
  hStack = ROOT.THStack("hStackTopologies", ";;" + ("Fraction" if normalize else "Number") + " of Combos (RF-subtracted)" + (" [%]" if normalize else ""))
  # get overall distribution
  hist = hists[""]
  hist.SetName("Overall" if normalize else "Total")
  hist.SetLineColor(colors[""])
  # get bin content of found and missing histograms in the order of the topologies in the overall histogram
  bins = [{"topology" : hist.GetXaxis().GetBinLabel(i)} for i in range(1, hist.GetNbinsX() + 1)]
  for bin in bins:
    for case in ["Found", "Missing"]:
      hist = hists[case]
      bin[case] = hist.GetBinContent(hist.GetXaxis().FindBin(bin["topology"]))
  # overlay missing and found histograms
  for case in ["Found", "Missing"]:
    hist = hists[case] = hists[""].Clone(case)
    for index, bin in enumerate(bins):
      hist.SetBinContent(index + 1, bin[case])
      assert bin["topology"] == hist.GetXaxis().GetBinLabel(index + 1)
    hist.SetLineColor(colors[case])
    hStack.Add(hist)
  hStack.Add(hists[""])  # draw on top
  # plot overall distribution
  canv = ROOT.TCanvas("topologies" + ("_norm" if normalize else ""))
  hStack.Draw("NOSTACK HIST")
  hStack.GetXaxis().SetRangeUser(0, maxNmbTopologies)
  hStack.SetMinimum(0)
  # add legend
  legend = canv.BuildLegend()
  # add labels that show number or fraction outside of plot range
  legend.AddEntry(ROOT.MakeNullPointer(ROOT.TObject), "Other topologies:", "")
  for case in cases:
    integralOtherTopos = hists[case].Integral(maxNmbTopologies, hists[case].GetNbinsX())
    legendEntry = legend.AddEntry(ROOT.MakeNullPointer(ROOT.TObject), "    " + str(round(integralOtherTopos)) + ("%" if normalize else " Combos"), "")
    legendEntry.SetTextColor(colors[case])
  canv.SaveAs(".pdf")
  return [bins[i]["topology"] for i in range(maxNmbTopologies)]  # return maxNmbTopologies largest topologies in overall distribution


def plotMissingMassSquared(topologies):
  inFileName = "pippippimpimpmiss.root"
  histBaseName = "MissingMassSquared/MissingMassSquared"
  rebinFactor = 100
  colorOffset = 1

  # get histograms
  inFile = ROOT.TFile(inFileName)
  cases = ["Found", "Missing", ""]
  hists = {topology : {case : inFile.Get(histBaseName + ("_" + case if case != "" else "") + "__" + topology) for case in cases} for topology in topologies}

  # overlay distributions for topologies
  hStacks = {case : ROOT.THStack("hStackMissingMassSquared" + case, ("Overall" if case == "" else case) + ";;Number of Combos (RF-subtracted)") for case in cases}
  for case in cases:
    for i, topology in enumerate(topologies):
      hist = hists[topology][case]
      hist.SetName(topology)
      hist.SetLineColor(i + colorOffset)
      hist.Rebin(rebinFactor)
      hStacks[case].Add(hist)
    canv = ROOT.TCanvas("mm2_topologies" + ("_" + case if case != "" else ""))
    hStacks[case].Draw("NOSTACK HIST")
    hStacks[case].GetXaxis().SetTitle(hists[topologies[0]][case].GetXaxis().GetTitle())
    hStacks[case].GetXaxis().SetTitleOffset(0)
    print(hStacks[case].GetXaxis().GetTitle())
    # add legend
    canv.BuildLegend()
    canv.Modified()
    canv.Update()
    canv.SaveAs(".pdf")

if __name__ == "__main__":
  topologies = plotTopologies(normalize = False)
  plotTopologies(normalize = True)
  plotMissingMassSquared(topologies)