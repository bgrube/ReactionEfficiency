#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import pandas
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetHistMinimumZero(True)


def plotTopologies(maxNmbTopologies = 10, normalize = False):
  inFileName = "pippippimpimpmiss_bggen_2017_01-ver03.root"
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
    print(f"{case} signal: {hist.GetBinContent(1)}{'%' if normalize else ' combos'}")
  hStack.Add(hists[""])  # draw on top
  print(f"Overall signal: {hists[''].GetBinContent(1)}{'%' if normalize else ' combos'}")
  # plot overall distribution
  canv = ROOT.TCanvas("justin_Proton_4pi_topologies" + ("_norm" if normalize else ""))
  hStack.Draw("NOSTACK HIST")
  hStack.GetXaxis().SetRangeUser(0, maxNmbTopologies)
  hStack.SetMinimum(0)
  # add legend
  legend = canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
  # add labels that show number or fraction outside of plot range
  legend.AddEntry(ROOT.MakeNullPointer(ROOT.TObject), "Other topologies:", "")
  for case in cases:
    integralOtherTopos = hists[case].Integral(maxNmbTopologies, hists[case].GetNbinsX())
    legendEntry = legend.AddEntry(ROOT.MakeNullPointer(ROOT.TObject), "    " + str(round(integralOtherTopos)) + ("%" if normalize else " Combos"), "")
    legendEntry.SetTextColor(colors[case])
  canv.SaveAs(".pdf")
  return [bins[i]["topology"] for i in range(maxNmbTopologies)]  # return maxNmbTopologies largest topologies in overall distribution


def plotMissingMassSquared(topologies):
  inFileName = "pippippimpimpmiss_bggen_2017_01-ver03.root"
  histBaseName = "MissingMassSquared/MissingMassSquared"
  rebinFactor = 100
  colorOffset = 1

  # get histograms
  inFile = ROOT.TFile(inFileName)
  cases = ["Found", "Missing", ""]
  hists = {topology : {case : inFile.Get(histBaseName + ("_" + case if case != "" else "") + "__" + topology) for case in cases} for topology in topologies}
  histsTotal = {case : inFile.Get(histBaseName + ("_" + case if case != "" else "")) for case in cases}

  # overlay distributions for topologies
  hStacks = {case : ROOT.THStack("hStackMissingMassSquaredTopologies" + case, ("Overall" if case == "" else case) + f";{hists[topologies[0]][case].GetXaxis().GetTitle()};Number of Combos (RF-subtracted)") for case in cases}
  for case in cases:
    # total distribution
    hist = histsTotal[case]
    hist.SetName("Total")
    hist.Rebin(rebinFactor)
    hist.SetFillColor(ROOT.kGray)
    hist.SetLineColor(ROOT.kGray)
    hStacks[case].Add(hist)
    # distribution for topologies
    for i, topology in enumerate(topologies):
      hist = hists[topology][case]
      hist.SetName(topology)
      hist.Rebin(rebinFactor)
      hist.SetLineColor(i + colorOffset)
      hStacks[case].Add(hist)
    canv = ROOT.TCanvas("justin_Proton_4pi_mm2_bggen_topologies" + ("_" + case if case != "" else ""))
    hStacks[case].Draw("NOSTACK HIST")
    # add legend
    canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
    canv.SaveAs(".pdf")


def overlayMissingMassSquared():
  inFileNames = ("pippippimpimpmiss.30370.root", "pippippimpimpmiss_bggen_2017_01-ver03.root")
  labels = ("Real data (scaled)", "bggen MC")
  histBaseName = "MissingMassSquared/MissingMassSquared"
  rebinFactor = 100

  # get histograms
  inFiles = [ROOT.TFile(inFileName) for inFileName in inFileNames]
  cases = ["Found", "Missing", ""]
  hists = [{case : inFile.Get(histBaseName + ("_" + case if case != "" else "")) for case in cases} for inFile in inFiles]

  # overlay real-data and bggen MC distributions
  hStacks = {case : ROOT.THStack("hStackMissingMassSquaredOverlay" + case, ("Total" if case == "" else case) + f";{hists[0][case].GetXaxis().GetTitle()};Number of Combos (RF-subtracted)") for case in cases}
  for case in cases:
    # normalize real data
    hist = hists[0][case]
    hist.Scale(hists[1][case].Integral() / hist.Integral())
    # set style
    hist.SetLineColor(ROOT.kGray)
    hist.SetFillColor(ROOT.kGray)
    hists[1][case].SetLineColor(ROOT.kRed + 1)
    for i, hist in enumerate(hists):
      hist = hists[i][case]
      hist.SetName(labels[i])
      hist.Rebin(rebinFactor)
      hStacks[case].Add(hist)
    # draw distributions
    canv = ROOT.TCanvas("justin_Proton_4pi_mm2_bggen_overlay" + ("_" + case if case != "" else ""))
    hStacks[case].Draw("NOSTACK HIST")
    # add legend
    canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
    canv.SaveAs(".pdf")


def plotMcTruthComparison():
  inFileName = "pippippimpimpmiss_bggen_2017_01-ver03.root"
  histBaseNames = (
    "MissingMassSquared/TruthDeltaPOverP",
    "MissingMassSquared/TruthDeltaTheta",
    "MissingMassSquared/TruthDeltaPhi")
  sigTopology = "2#pi^{#plus}2#pi^{#minus}p"
  bkgTopology = "2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0}]"
  colors = {
    "Found"   : ROOT.kGreen + 2,
    "Missing" : ROOT.kRed + 1,
    ""        : ROOT.kGray
  }

  # get histograms
  inFile = ROOT.TFile(inFileName)
  cases = ["", "Found", "Missing"]
  histsTot = [{case : inFile.Get(histBaseName + ("_" + case if case != "" else "")) for case in cases} for histBaseName in histBaseNames]
  histsSig = [{case : inFile.Get(histBaseName + ("_" + case if case != "" else "") + "__" + sigTopology) for case in cases} for histBaseName in histBaseNames]
  histsBkg = [{case : inFile.Get(histBaseName + ("_" + case if case != "" else "") + "__" + bkgTopology) for case in cases} for histBaseName in histBaseNames]
  histsBkgTot = []
  for i, hists in enumerate(histsTot):
    histsBkgTot.append({})
    for case, histTot in hists.items():
      histSig = histsSig[i][case]
      histBkg = histTot.Clone(histSig.GetName().split("__")[0] + "__bkg")
      histBkg.Add(histSig, -1)
      histsBkgTot[i][case] = histBkg

  # draw histograms
  for hists in histsTot + histsSig + histsBkg + histsBkgTot:
    distrName = hists[""].GetName()
    hStack = ROOT.THStack("hStack" + distrName, f";{hists[''].GetXaxis().GetTitle()};Number of Combos (RF-subtracted)")
    canv = ROOT.TCanvas("justin_Proton_4pi_mm2_bggen_mctruthcomp_" + distrName, "", 600, 600)
    for case, hist in hists.items():
      hist.SetName(case if case != "" else "Total")
      hist.SetTitle("")
      hist.SetLineColor(colors[case])
      if case == "":
        hist.SetFillColor(colors[case])
      hStack.Add(hist)
    hStack.Draw("NOSTACK HIST")
    hStack.GetYaxis().SetTitleOffset(1.5)
    if "TruthDeltaTheta" in distrName:
      hStack.GetXaxis().SetRangeUser(-100, 100)
    # add legend
    canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
    canv.SaveAs(".pdf")


def overlayTopologies(
  inFileName,
  histName,
  maxNmbTopologies = 10,
  xRange           = None):
  # get histograms
  inFile = ROOT.TFile(inFileName)
  hist3D = inFile.Get(histName)
  xAxisTitle = hist3D.GetXaxis().GetTitle()
  caseAxis = hist3D.GetYaxis()
  cases = [caseAxis.GetBinLabel(i + 1) for i in range(0, caseAxis.GetNbins())]
  # topoAxis = hist3D.GetZaxis()
  # topologies = [topoAxis.GetBinLabel(i + 1) for i in range(0, topoAxis.GetNbins())]
  # assert topologies[0] == "Total"
  for case in cases:
    hStack = ROOT.THStack(f"{hist3D.GetName()}_{case}", f"{case};{xAxisTitle};Number of Combos (RF-subtracted)")
    caseBin = caseAxis.FindBin(case)
    caseAxis.SetRange(caseBin, caseBin)
    hist2D = hist3D.Project3D("ZX_" + case)
    hist2D.LabelsDeflate("Y")
    hist2D.LabelsOption(">", "Y")
    topoAxis = hist2D.GetYaxis()
    topologies = [topoAxis.GetBinLabel(i + 1) for i in range(0, topoAxis.GetNbins())]
    # topoHistSorted = []
    # for yBin in range(1, topoAxis.GetNbins() + 1):
    #   nmbCombos = 0
    #   for xBin in range(1, hist2D.GetXaxis().GetNbins() + 1):
    #     nmbCombos += hist2D.GetBinContent(xBin, yBin)
    #   topoHistSorted.append((topoAxis.GetBinLabel(yBin), nmbCombos))
    # print(case, topoHistSorted[:maxNmbTopologies + 1])
    assert topologies[0] == "Total"
    hists = []
    # overlay distributions for topologies
    for i in range(0, maxNmbTopologies + 1):
      hist1D = hist2D.ProjectionX(f"{hist3D.GetName()}_{case}_{topologies[i]}", i + 1, i + 1)
      hist1D.SetTitle(topologies[i])
      if i == 0:
        hist1D.SetLineColor(ROOT.kGray)
        hist1D.SetFillColor(ROOT.kGray)
      else:
        hist1D.SetLineColor(i)
      hists.append(hist1D)
      hStack.Add(hist1D)
    # draw distributions
    canv = ROOT.TCanvas(f"justin_Proton_4pi_{hist3D.GetName()}_bggen_topologies_{case}")
    hStack.Draw("NOSTACK HIST")
    xAxis = hStack.GetXaxis()
    if not xRange is None:
      xAxis.SetRangeUser(*xRange)
    # add legend
    canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
    if hStack.GetMinimum() < 0 and hStack.GetMaximum() > 0:
      line = ROOT.TLine()
      line.SetLineStyle(ROOT.kDashed)
      line.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 0, xAxis.GetBinUpEdge(xAxis.GetLast()), 0)
    canv.SaveAs(".pdf")


def overlayTopologies2(
  inFileName,
  treeName,
  branchName,
  xTitle,
  binning,
  maxNmbTopologies = 10):
  inputData = ROOT.RDataFrame(treeName, inFileName)
  # get topologies with largest number of combos for total data set
  # print(df.GetColumnType("ThrownTopology"))
  toposToPlot = ["Total"]
  topos = [str(topo) for topo in inputData.AsNumpy(["ThrownTopology"])["ThrownTopology"]]
  if len(topos) > 1:
    weights = inputData.AsNumpy(["AccidWeightFactor"])["AccidWeightFactor"]
    pdf = pandas.DataFrame({"ThrownTopology" : topos, "AccidWeightFactor" : weights})
    topoHistSorted = pdf.groupby("ThrownTopology")["AccidWeightFactor"].sum().sort_values(ascending = False)
    toposToPlot += list(topoHistSorted[:maxNmbTopologies].index)
  # print(toposToPlot)
  filterCases = {
    "Total"   : "true",
    "Found"   : "TrackFound == true",
    "Missing" : "TrackFound == false"
  }
  for case, caseFilter in filterCases.items():
    caseData = inputData.Filter(caseFilter)
    # get topologies with largest number of combos for given case
    toposToPlot = ["Total"]
    topos = [str(topo) for topo in caseData.AsNumpy(["ThrownTopology"])["ThrownTopology"]]
    if len(topos) > 1:
      weights = caseData.AsNumpy(["AccidWeightFactor"])["AccidWeightFactor"]
      pdf = pandas.DataFrame({"ThrownTopology" : topos, "AccidWeightFactor" : weights})
      topoHistSorted = pdf.groupby("ThrownTopology")["AccidWeightFactor"].sum().sort_values(ascending = False)
      # print(case, topoHistSorted[:maxNmbTopologies])
      toposToPlot += list(topoHistSorted[:maxNmbTopologies].index)
    hStack = ROOT.THStack(f"{branchName}_{case}", f"{case};{xTitle};Number of Combos (RF-subtracted)")
    hists = []
    # overlay distributions for topologies
    for index, topo in enumerate(toposToPlot):
      topoData = caseData.Filter(f'ThrownTopology.GetString() == "{topo}"' if topo != "Total" else "true")
      hist = topoData.Histo1D((f"{branchName}_{case}_{topo}", topo, *binning), branchName, "AccidWeightFactor")
      if topo == "Total":
        hist.SetLineColor(ROOT.kGray)
        hist.SetFillColor(ROOT.kGray)
      else:
        hist.SetLineColor(index)
      hists.append(hist)
      hStack.Add(hist.GetPtr())
    # draw distributions
    canv = ROOT.TCanvas(f"justin_Proton_4pi_{branchName}_bggen_topologies_{case}_new")
    hStack.Draw("NOSTACK HIST")
    xAxis = hStack.GetXaxis()
    # add legend
    canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
    # draw zero line, if necessary
    if hStack.GetMinimum() < 0 and hStack.GetMaximum() > 0:
      line = ROOT.TLine()
      line.SetLineStyle(ROOT.kDashed)
      line.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 0, xAxis.GetBinUpEdge(xAxis.GetLast()), 0)
    canv.SaveAs(".pdf")


if __name__ == "__main__":
  ROOT.gROOT.LoadMacro("~/rootlogon.C")
  ROOT.gROOT.ForceStyle()
  ROOT.gStyle.SetCanvasDefW(600)
  ROOT.gStyle.SetCanvasDefH(400)
  ROOT.gStyle.SetPalette(ROOT.kBird)
  ROOT.gStyle.SetLegendFillColor(ROOT.kWhite)
  ROOT.gStyle.SetLegendBorderSize(1)
  # ROOT.gStyle.SetOptStat("ni")  # show only name and integral
  ROOT.gStyle.SetOptStat("i")  # show only integral
  ROOT.gStyle.SetStatFormat("8.8g")
  ROOT.gStyle.SetTitleColor(1, "X")  # fix that for some mysterious reason x-axis titles of 2D plots and graphs are white

  inFileName = "pippippimpimpmiss.root"
  # topologies = plotTopologies(normalize = False)
  # plotTopologies(normalize = True)
  # plotMissingMassSquared(topologies)
  # overlayMissingMassSquared()
  # plotMcTruthComparison()
  overlayTopologies(inFileName, "MissingMassSquared/NmbUnusedShowers",    xRange = (-0.5, 10.5))
  overlayTopologies(inFileName, "MissingMassSquared/EnergyUnusedShowers", xRange = (0,    6))

  inFileName = "pippippimpimpmiss_flatTree.root"
  treeName = "pippippimpimpmiss"
  overlayTopologies2(inFileName, treeName, "NmbUnusedShowers",    "Number of Unused Showers",   binning = (11, -0.5, 10.5))
  overlayTopologies2(inFileName, treeName, "EnergyUnusedShowers", "Unused Shower Energy (GeV)", binning = (60,  0,    6))
