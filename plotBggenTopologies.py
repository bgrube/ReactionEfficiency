#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


from collections.abc import Iterable
import pandas
import ROOT

ROOT.gROOT.SetBatch(True)


FILTER_CASES = {
  "Total"   : "true",
  "Found"   : "TrackFound == true",
  "Missing" : "TrackFound == false"
}


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


def drawHistogram(
  inFileName,
  histName,
  rebinFactor = 100,
  pdfFileNamePrefix = "justin_Proton_4pi_",
  pdfFileNameSuffix = ""
):
  # get histogram
  inFile = ROOT.TFile(inFileName)
  hist = inFile.Get(histName)
  hist.Rebin(rebinFactor)
  # draw histogram
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{hist.GetName()}{pdfFileNameSuffix}")
  hist.Draw("HIST")
  canv.SaveAs(".pdf")


def getHistND(
  inputData,       # RDataFrame
  variables,       # tuple of either column names or tuples with new column definitions; defines dimension of histogram
  binning,
  titles,          # tuple with axis titles
  weightVariable,  # may be None (= no weighting), string with column name, or tuple with new column definition
  filterExpression
):
  histDim = len(variables)
  assert 1 <= histDim <= 2, "currently, only 1D and 2D histograms are supported"
  if filterExpression:
    # apply additional filters
    inputData = inputData.Filter(filterExpression)
  columnNames = [None,] * len(variables)
  for index, variable in enumerate(variables):
    if isinstance(variable, str):
      # use existing variable column
      columnNames[index] = variable
    elif isinstance(variable, Iterable):
      # create new variable column
      inputData  = inputData.Define(variable[0], variable[1])
      columnNames[index] = variable[0]
    assert columnNames[index] is not None, f"failed to get column name for variable '{variable}'"
  if not isinstance(weightVariable, str) and isinstance(weightVariable, Iterable):
    # create new weight column
    inputData = inputData.Define(weightVariable[0], weightVariable[1])
  # create histogram
  hist = None
  histDef = ("_vs_".join(columnNames), ";" + ";".join(titles), *binning)
  # get member function to create histogram
  HistoND = getattr(inputData, "Histo1D") if histDim == 1 else getattr(inputData, "Histo2D")
  if weightVariable is None:
    hist = HistoND(histDef, *columnNames)
  elif isinstance(weightVariable, str):
    # use existing weight column
    hist = HistoND(histDef, *columnNames, weightVariable)
  elif isinstance(weightVariable, Iterable):
    # use new weight column
    hist = HistoND(histDef, *columnNames, weightVariable[0])
  assert hist is not None, f"failed to create histogram for weight variable '{weightVariable}'"
  return hist


# plot distribution of variable defined by `variable`
def plotVariable(
  inputData,  # RDataFrame
  variable,  # may be column name, or tuple with new column definition
  binning,
  xTitle,
  yTitle = "Number of Combos (RF-subtracted)",
  weightVariable = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
  pdfFileNamePrefix = "justin_Proton_4pi_",
  pdfFileNameSuffix = "",
  additionalFilter  = None
):
  hist = getHistND(inputData, (variable,), binning, (xTitle, yTitle), weightVariable, additionalFilter)
  # draw distributions
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{hist.GetName()}{pdfFileNameSuffix}")
  hist.Draw("HIST")
  hist.GetYaxis().SetTitleOffset(1.5)
  # draw zero line, if necessary
  xAxis = hist.GetXaxis()
  if hist.GetMinimum() < 0 and hist.GetMaximum() > 0:
    line = ROOT.TLine()
    line.SetLineStyle(ROOT.kDashed)
    line.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 0, xAxis.GetBinUpEdge(xAxis.GetLast()), 0)
  canv.SaveAs(".pdf")


# overlay distributions of variable defined by `variable` for Total, Found, and Missing cases
def overlayCases(
  inFileName,
  treeName,
  variable,
  xTitle,
  binning,
  additionalFilter  = None,
  pdfFileNamePrefix = "justin_Proton_4pi_",
  pdfFileNameSuffix = ""
):
  inputData = ROOT.RDataFrame(treeName, inFileName) if additionalFilter is None else ROOT.RDataFrame(treeName, inFileName).Filter(additionalFilter)
  colorCases = {
    "Total"   : ROOT.kGray,
    "Found"   : ROOT.kGreen + 2,
    "Missing" : ROOT.kRed + 1
  }
  # overlay distributions for cases
  hStack = ROOT.THStack(f"{variable}", f";{xTitle};Number of Combos (RF-subtracted)")
  hists = []
  for case in FILTER_CASES.keys():
    caseData = inputData.Filter(FILTER_CASES[case])
    hist = caseData.Histo1D((f"{variable}_{case}", case, *binning), variable, "AccidWeightFactor")
    hist.SetLineColor(colorCases[case])
    if case == "Total":
      hist.SetFillColor(colorCases[case])
    hists.append(hist)
    hStack.Add(hist.GetPtr())
  # draw distributions
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{variable}_cases{pdfFileNameSuffix}")
  hStack.Draw("NOSTACK HIST")
  hStack.GetYaxis().SetTitleOffset(1.5)
  # add legend
  canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
  # draw zero line, if necessary
  xAxis = hStack.GetXaxis()
  if hStack.GetMinimum() < 0 and hStack.GetMaximum() > 0:
    line = ROOT.TLine()
    line.SetLineStyle(ROOT.kDashed)
    line.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 0, xAxis.GetBinUpEdge(xAxis.GetLast()), 0)
  canv.SaveAs(".pdf")


# plot distribution of variable defined by `variable` for overall data sample
# for bggen sample: overlay distributions for the `maxNmbTopologies` topologies with the largest number of combos
def overlayTopologies(
  inFileName,
  treeName,
  variable,
  xTitle,
  binning,
  additionalFilter  = None,
  pdfFileNamePrefix = "justin_Proton_4pi_",
  maxNmbTopologies  = 10
):
  inputData = ROOT.RDataFrame(treeName, inFileName) if additionalFilter is None else ROOT.RDataFrame(treeName, inFileName).Filter(additionalFilter)
  # # get topologies with largest number of combos for total data set
  # # print(df.GetColumnType("ThrownTopology"))
  # toposToPlot = ["Total"]
  # topos = [str(topo) for topo in inputData.AsNumpy(["ThrownTopology"])["ThrownTopology"]]
  # if len(topos) > 1:
  #   weights = inputData.AsNumpy(["AccidWeightFactor"])["AccidWeightFactor"]
  #   pdf = pandas.DataFrame({"ThrownTopology" : topos, "AccidWeightFactor" : weights})
  #   topoHistSorted = pdf.groupby("ThrownTopology")["AccidWeightFactor"].sum().sort_values(ascending = False)
  #   toposToPlot += list(topoHistSorted[:maxNmbTopologies].index)
  # # print(toposToPlot)
  for case in FILTER_CASES.keys():
    caseData = inputData.Filter(FILTER_CASES[case])
    # get topologies with largest number of combos for given case
    toposToPlot = ["Total"]
    topos = [str(topo) for topo in caseData.AsNumpy(["ThrownTopology"])["ThrownTopology"]]
    if len(topos) > 1:
      weights = caseData.AsNumpy(["AccidWeightFactor"])["AccidWeightFactor"]
      pdf = pandas.DataFrame({"ThrownTopology" : topos, "AccidWeightFactor" : weights})
      topoHistSorted = pdf.groupby("ThrownTopology")["AccidWeightFactor"].sum().sort_values(ascending = False)
      # print(case, topoHistSorted[:maxNmbTopologies])
      toposToPlot += list(topoHistSorted[:maxNmbTopologies].index)
    hStack = ROOT.THStack(f"{variable}_{case}", f"{case};{xTitle};Number of Combos (RF-subtracted)")
    hists = []
    # overlay distributions for topologies
    for index, topo in enumerate(toposToPlot):
      topoData = caseData.Filter(f'ThrownTopology.GetString() == "{topo}"' if topo != "Total" else "true")
      hist = topoData.Histo1D((f"{variable}_{case}_{topo}", topo, *binning), variable, "AccidWeightFactor")
      if topo == "Total":
        hist.SetLineColor(ROOT.kGray)
        hist.SetFillColor(ROOT.kGray)
      else:
        hist.SetLineColor(index)
      hists.append(hist)
      hStack.Add(hist.GetPtr())
    # draw distributions
    canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{variable}_bggen_topologies_{case}")  #TODO replace "_bggen_topologies_" with configurable string
    hStack.Draw("NOSTACK HIST")
    hStack.GetYaxis().SetTitleOffset(1.5)
    # add legend
    canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
    # draw zero line, if necessary
    xAxis = hStack.GetXaxis()
    if hStack.GetMinimum() < 0 and hStack.GetMaximum() > 0:
      line = ROOT.TLine()
      line.SetLineStyle(ROOT.kDashed)
      line.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 0, xAxis.GetBinUpEdge(xAxis.GetLast()), 0)
    canv.SaveAs(".pdf")


if __name__ == "__main__":
  ROOT.gROOT.LoadMacro("~/rootlogon.C")
  ROOT.gROOT.ForceStyle()
  ROOT.gStyle.SetCanvasDefW(600)
  ROOT.gStyle.SetCanvasDefH(600)
  ROOT.gStyle.SetPalette(ROOT.kBird)
  ROOT.gStyle.SetLegendFillColor(ROOT.kWhite)
  ROOT.gStyle.SetLegendBorderSize(1)
  # ROOT.gStyle.SetOptStat("ni")  # show only name and integral
  ROOT.gStyle.SetOptStat("i")  # show only integral
  ROOT.gStyle.SetStatFormat("8.8g")
  ROOT.gStyle.SetTitleColor(1, "X")  # fix that for some mysterious reason x-axis titles of 2D plots and graphs are white

  histFileName = "pippippimpimpmiss.root"
  # topologies = plotTopologies(normalize = False)
  # plotTopologies(normalize = True)
  # overlayMissingMassSquared()

  treeFileName = "pippippimpimpmiss_flatTree.root"
  treeName = "pippippimpimpmiss"

  overlayTopologies(treeFileName, treeName, "NmbUnusedShowers",            xTitle = "Number of Unused Showers",             binning = (11, -0.5, 10.5))
  overlayTopologies(treeFileName, treeName, "EnergyUnusedShowers",         xTitle = "Unused Shower Energy (GeV)",           binning = (60, 0, 6))
  overlayTopologies(treeFileName, treeName, "MissingMassSquared_Measured", xTitle = "Missing Mass Squared (GeV/c^{2})^{2}", binning = (50, -0.5, 4.5))

  filterTopologies = {
    ""                                             : None,
    "__2#pi^{#plus}2#pi^{#minus}p"                 : 'ThrownTopology.GetString() == "2#pi^{#plus}2#pi^{#minus}p"',
    "__2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0}]" : 'ThrownTopology.GetString() == "2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0}]"',
    "__bkg"                                        : 'ThrownTopology.GetString() != "2#pi^{#plus}2#pi^{#minus}p"'
  }
  for suffix, filter in filterTopologies.items():
    overlayCases(treeFileName, treeName, "TruthDeltaP",      xTitle = "#it{p}^{miss}_{truth} #minus #it{p}^{miss}_{kin. fit} (GeV/c)",                      binning = (600, -6, 6),     additionalFilter = filter, pdfFileNameSuffix = suffix)
    overlayCases(treeFileName, treeName, "TruthDeltaPOverP", xTitle = "(#it{p}^{miss}_{truth} #minus #it{p}^{miss}_{kin. fit}) / #it{p}^{miss}_{kin. fit}", binning = (500, -2, 2),     additionalFilter = filter, pdfFileNameSuffix = suffix)
    overlayCases(treeFileName, treeName, "TruthDeltaTheta",  xTitle = "#it{#theta}^{miss}_{truth} #minus #it{#theta}^{miss}_{kin. fit} (deg)",              binning = (200, -100, 100), additionalFilter = filter, pdfFileNameSuffix = suffix)
    overlayCases(treeFileName, treeName, "TruthDeltaPhi",    xTitle = "#it{#phi}^{miss}_{truth} #minus #it{#phi}^{miss}_{kin. fit} (deg)",                  binning = (360, -180, 180), additionalFilter = filter, pdfFileNameSuffix = suffix)

  inputData = ROOT.RDataFrame(treeName, treeFileName)
  sideBandArgs = {"yTitle" : "Number of Combos (RF-Sideband)", "weightVariable" : ("AccidWeightFactorSb", "1 - AccidWeightFactor"), "pdfFileNameSuffix" : "_Sb"}
  plotVariable(inputData, ("MissingMass_Measured", "sqrt(MissingMassSquared_Measured)"), xTitle = "Missing Mass (GeV/c^{2})",             binning = (100, 0, 2))
  plotVariable(inputData, ("MissingMass_Measured", "sqrt(MissingMassSquared_Measured)"), xTitle = "Missing Mass (GeV/c^{2})",             binning = (100, 0, 2), **sideBandArgs)
  plotVariable(inputData, "MissingMassSquared_Measured",                                 xTitle = "Missing Mass Squared (GeV/c^{2})^{2}", binning = (225, -0.5, 4))
  plotVariable(inputData, "MissingMassSquared_Measured",                                 xTitle = "Missing Mass Squared (GeV/c^{2})^{2}", binning = (225, -0.5, 4), **sideBandArgs)
  plotVariable(inputData, "AccidWeightFactor",                                           xTitle = "RF Weight",                            binning = (1000, -2, 2), weightVariable = None)
