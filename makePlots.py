#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


from collections.abc import Iterable
import pandas
import ROOT


# simplified version of Paul's criteria
TRACK_FOUND_FILTER = "(" \
  + "(NmbUnusedTracks == 1)" \
  + " and ((MissingProtonTheta < 5) or (abs(UnusedDeltaPhi[0]) <= 30))" \
  + " and (abs(UnusedDeltaTheta[0]) <= 30)" \
  + " and (abs(UnusedDeltaPOverP[0]) <= 0.6)" \
  + ")"

# filter expressions for track-found cases
FILTER_CASES = {
  "Total"   : "(true)",
  "Found"   : "(TrackFound == true)",
  "Missing" : "(TrackFound == false)"
}


def setupPlotStyle():
  #TODO remove dependency from external file or add file to repo
  ROOT.gROOT.LoadMacro("~/rootlogon.C")
  ROOT.gROOT.ForceStyle()
  ROOT.gStyle.SetCanvasDefW(600)
  ROOT.gStyle.SetCanvasDefH(600)
  ROOT.gStyle.SetPalette(ROOT.kBird)
  # ROOT.gStyle.SetPalette(ROOT.kViridis)
  ROOT.gStyle.SetLegendFillColor(ROOT.kWhite)
  ROOT.gStyle.SetLegendBorderSize(1)
  # ROOT.gStyle.SetOptStat("ni")  # show only name and integral
  ROOT.gStyle.SetOptStat("i")  # show only integral
  ROOT.gStyle.SetStatFormat("8.8g")
  ROOT.gStyle.SetTitleColor(1, "X")  # fix that for some mysterious reason x-axis titles of 2D plots and graphs are white
  ROOT.gStyle.SetTitleOffset(1.35, "Y")


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


# simple generic function to plot a histogram in a ROOT file
def drawHistogram(
  inFileName,
  histName,
  rebinFactor = 1,  # if integer -> rebin x-axis, if iterable of integers rebin axes
  drawOption = "HIST",
  pdfFileNamePrefix = "justin_Proton_4pi_",
  pdfFileNameSuffix = ""
):
  # get histogram
  inFile = ROOT.TFile(inFileName)
  hist = inFile.Get(histName)
  if isinstance(hist, ROOT.TH2):
    if isinstance(rebinFactor, int):
      hist.RebinX(rebinFactor)
    elif isinstance(rebinFactor, Iterable):
      hist.Rebin2D(rebinFactor[0], rebinFactor[1])
  elif isinstance(hist, ROOT.TH1):
    hist.Rebin(rebinFactor)
  # draw histogram
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{hist.GetName()}{pdfFileNameSuffix}")
  hist.Draw(drawOption)
  canv.SaveAs(".pdf")


def getHistND(
  inputData,   # RDataFrame
  variables,   # variable(s) to plot; is tuple of either column names or tuples with new column definitions; defines dimension of histogram
  axisTitles,  # semicolon-separated list
  binning,     # tuple with binning definition
  weightVariable   = None,  # may be None (= no weighting), string with column name, or tuple with new column definition
  filterExpression = None,
  histNameSuffix   = "",
  histTitle        = ""
):
  histDim = len(variables)
  assert 1 <= histDim <= 2, "currently, only 1D and 2D histograms are supported"
  # apply additional filters, if defined
  data = inputData.Filter(filterExpression) if filterExpression else inputData
  columnNames = [None,] * len(variables)
  for index, variable in enumerate(variables):
    if isinstance(variable, str):
      # use existing variable column
      columnNames[index] = variable
    elif isinstance(variable, Iterable):
      # create new variable column
      data  = data.Define(variable[0], variable[1])
      columnNames[index] = variable[0]
    assert columnNames[index] is not None, f"failed to get column name for variable '{variable}'"
  if not isinstance(weightVariable, str) and isinstance(weightVariable, Iterable):
    # create new weight column
    data = data.Define(weightVariable[0], weightVariable[1])
  # create histogram
  hist = None
  histDef = ("_".join(columnNames + ([histNameSuffix] if histNameSuffix else [])), f"{histTitle};{axisTitles}", *binning)  # histogram definition
  # get member function to create histogram
  HistoND = getattr(data, "Histo1D") if histDim == 1 else getattr(data, "Histo2D")
  if not weightVariable:
    hist = HistoND(histDef, *columnNames)
  elif isinstance(weightVariable, str):
    # use existing weight column
    hist = HistoND(histDef, *columnNames, weightVariable)
  elif isinstance(weightVariable, Iterable):
    # use new weight column
    hist = HistoND(histDef, *columnNames, weightVariable[0])
  assert hist is not None, f"failed to create histogram for weight variable '{weightVariable}'"
  return hist


def setDefaultYAxisTitle(
  axisTitles,  # semicolon-separated list
  defaultYTitle = "Number of Combos (RF-subtracted)"
):
  titles = axisTitles.split(";")
  if (len(titles) == 1):
    return titles[0] + ";" + defaultYTitle
  elif (len(titles) == 0):
    return ";" + defaultYTitle
  else:
    return axisTitles


# plot 1D distribution of given variable
def plot1D(
  inputData,   # RDataFrame
  variable,    # variable to plot; may be column name, or tuple with new column definition
  axisTitles,  # semicolon-separated list
  binning,     # tuple with binning definition
  weightVariable = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
  pdfFileNamePrefix = "justin_Proton_4pi_",
  pdfFileNameSuffix = "",
  additionalFilter  = None
):
  hist = getHistND(inputData, (variable,), setDefaultYAxisTitle(axisTitles), binning, weightVariable, additionalFilter)
  # draw distributions
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{hist.GetName()}{pdfFileNameSuffix}")
  hist.Draw("HIST")
  # draw zero line, if necessary
  xAxis = hist.GetXaxis()
  if hist.GetMinimum() < 0 and hist.GetMaximum() > 0:
    line = ROOT.TLine()
    line.SetLineStyle(ROOT.kDashed)
    line.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 0, xAxis.GetBinUpEdge(xAxis.GetLast()), 0)
  canv.SaveAs(".pdf")


# plot 2D distribution of given variables
def plot2D(
  inputData,   # RDataFrame
  xVariable,   # x variable to plot; may be column name, or tuple with new column definition
  yVariable,   # y variable to plot; may be column name, or tuple with new column definition
  axisTitles,  # semicolon-separated list
  binning,     # tuple with binning definition
  weightVariable = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
  pdfFileNamePrefix = "justin_Proton_4pi_",
  pdfFileNameSuffix = "",
  additionalFilter  = None
):
  hist = getHistND(inputData, (xVariable, yVariable), axisTitles, binning, weightVariable, additionalFilter)
  # draw distributions
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{hist.GetName()}{pdfFileNameSuffix}")
  hist.Draw("COLZ")
  canv.SaveAs(".pdf")


# overlay distributions of variable defined by `variable` for Total, Found, and Missing cases
def overlayCases(
  inputData,   # RDataFrame
  variable,    # variable to plot; may be column name, or tuple with new column definition
  axisTitles,  # semicolon-separated list
  binning,     # tuple with binning definition
  additionalFilter  = None,
  pdfFileNamePrefix = "justin_Proton_4pi_",
  pdfFileNameSuffix = ""
):
  data = inputData.Filter(additionalFilter) if additionalFilter else inputData
  colorCases = {
    "Total"   : ROOT.kGray,
    "Found"   : ROOT.kGreen + 2,
    "Missing" : ROOT.kRed + 1
  }
  # overlay distributions for cases
  hStack = ROOT.THStack(f"{variable}", ";" + setDefaultYAxisTitle(axisTitles))
  hists = []
  for case in FILTER_CASES.keys():
    hist = getHistND(data, (variable,), setDefaultYAxisTitle(axisTitles), binning, "AccidWeightFactor", FILTER_CASES[case],
                     histNameSuffix = case, histTitle = case)
    hist.SetLineColor(colorCases[case])
    if case == "Total":
      hist.SetFillColor(colorCases[case])
    hists.append(hist)
    hStack.Add(hist.GetPtr())
  # draw distributions
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{variable}_cases{pdfFileNameSuffix}")
  hStack.Draw("NOSTACK HIST")
  # add legend
  canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
  # draw zero line, if necessary
  xAxis = hStack.GetXaxis()
  if hStack.GetMinimum() < 0 and hStack.GetMaximum() > 0:
    line = ROOT.TLine()
    line.SetLineStyle(ROOT.kDashed)
    line.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 0, xAxis.GetBinUpEdge(xAxis.GetLast()), 0)
  canv.SaveAs(".pdf")


def getTopologyHist(
  inputData  # RDataFrame
):
  # print(inputData.GetColumnType("ThrownTopology"))
  topos = [str(topo) for topo in inputData.AsNumpy(["ThrownTopology"])["ThrownTopology"]]
  if len(topos) > 1:
    weights = inputData.AsNumpy(["AccidWeightFactor"])["AccidWeightFactor"]
    df = pandas.DataFrame({"ThrownTopology" : topos, "AccidWeightFactor" : weights})
    topoHistSorted = df.groupby("ThrownTopology")["AccidWeightFactor"].sum().sort_values(ascending = False)
    # print(type(topoHistSorted), topoHistSorted)
    # print("INDEX", type(topoHistSorted.index), topoHistSorted.index)
    # print("VAKUES", type(topoHistSorted.values), topoHistSorted.values)
    return (list(topoHistSorted.index), topoHistSorted.values)


def plotTopologyHist(
  inputData,  # RDataFrame
  normalize         = False,
  maxNmbTopologies  = 10,
  additionalFilter  = None,
  pdfFileNamePrefix = "justin_Proton_4pi_",
  pdfFileNameSuffix = ""
):
  data = inputData.Filter(additionalFilter) if additionalFilter else inputData
  colorCases = {
    # "Total"   : ROOT.kGray,
    # "Found"   : ROOT.kGreen + 2,
    # "Missing" : ROOT.kRed + 1
    "Total"   : ROOT.kBlue,
    "Found"   : ROOT.kGreen + 2,
    "Missing" : ROOT.kRed + 1
  }
  # get histogram data
  topoHists = {}
  for case in FILTER_CASES.keys():
    caseData = data.Filter(FILTER_CASES[case])
    topoHists[case] = getTopologyHist(caseData)
  # overlay distributions for cases
  hStack = ROOT.THStack(f"topologies",  ";;" + ("Fraction" if normalize else "Number") + " of Combos (RF-subtracted)" + (" [%]" if normalize else ""))
  hists = {}
  topoLabels = topoHists["Total"][0]
  for case in FILTER_CASES.keys():
    hist = ROOT.TH1F(f"topologies_{case}{'_norm' if normalize else ''}", case, 1, 0, 1)
    histValues = dict(zip(topoHists[case][0], topoHists[case][1]))
    # set bin content of histogram
    for binLabel in topoLabels:
      hist.Fill(binLabel, histValues[binLabel] if binLabel in histValues else 0)
    if normalize:
      hist.Scale(100 / hist.Integral())
    hist.SetLineColor(colorCases[case])
    # if case == "Total":
    #   hist.SetFillColor(colorCases[case])
    hists[case] = hist
    hStack.Add(hist)
    print(f"{case} signal: {hist.GetBinContent(1)}{'%' if normalize else ' combos'}")
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}topologies{'_norm' if normalize else ''}{pdfFileNameSuffix}")
  hStack.Draw("NOSTACK HIST")
  hStack.GetXaxis().SetRangeUser(0, maxNmbTopologies)
  # # add legend
  legend = canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
  # add labels that show number or fraction outside of plot range
  legend.AddEntry(ROOT.MakeNullPointer(ROOT.TObject), "Other topologies:", "")
  for case in FILTER_CASES.keys():
    integralOtherTopos = hists[case].Integral(maxNmbTopologies, hists[case].GetNbinsX())
    legendEntry = legend.AddEntry(ROOT.MakeNullPointer(ROOT.TObject), "    " + str(round(integralOtherTopos)) + ("%" if normalize else " Combos"), "")
    legendEntry.SetTextColor(colorCases[case])
  canv.SaveAs(".pdf")


# plot distribution of variable defined by `variable` for overall data sample
# for bggen sample: overlay distributions for the `maxNmbTopologies` topologies with the largest number of combos
def overlayTopologies(
  inputData,   # RDataFrame
  variable,    # variable to plot; may be column name, or tuple with new column definition
  axisTitles,  # semicolon-separated list
  binning,     # tuple with binning definition
  additionalFilter  = None,
  maxNmbTopologies  = 10,
  pdfFileNamePrefix = "justin_Proton_4pi_",
  pdfFileNameSuffix = "_bggen_topologies"
):
  data = inputData.Filter(additionalFilter) if additionalFilter else inputData
  # # get topologies with largest number of combos for total data set
  # toposToPlot, _ = getTopologyHist(data)
  # toposToPlot = ["Total"] + toposToPlot[:maxNmbTopologies]
  for case in FILTER_CASES.keys():
    caseData = data.Filter(FILTER_CASES[case])
    # get topologies with largest number of combos for given case
    toposToPlot, _ = getTopologyHist(caseData)
    toposToPlot = ["Total"] + toposToPlot[:maxNmbTopologies]
    hStack = ROOT.THStack(f"{variable}_{case}", f"{case};{setDefaultYAxisTitle(axisTitles)}")
    hists = []
    # overlay distributions for topologies
    for index, topo in enumerate(toposToPlot):
      hist = getHistND(caseData, (variable,), setDefaultYAxisTitle(axisTitles), binning, "AccidWeightFactor", (f'ThrownTopology.GetString() == "{topo}"' if topo != "Total" else "true"),
                       histNameSuffix = f"{case}_{topo}", histTitle = topo)
      if topo == "Total":
        hist.SetLineColor(ROOT.kGray)
        hist.SetFillColor(ROOT.kGray)
      else:
        hist.SetLineColor(index)
      hists.append(hist)
      hStack.Add(hist.GetPtr())
    # draw distributions
    canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{variable}{pdfFileNameSuffix}_{case}")
    hStack.Draw("NOSTACK HIST")
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
  ROOT.gROOT.SetBatch(True)
  ROOT.EnableImplicitMT(20)  # activate implicit multi-threading for RDataFrame; disable using ROOT.DisableImplicitMT()
  setupPlotStyle()

  # overlayMissingMassSquared()

  histFileName = "pippippimpimpmiss.root"
  treeFileName = "pippippimpimpmiss_flatTree.root"
  # histFileName = "pippippimpimpmiss_bggen_2017_01-ver03.root"
  # treeFileName = "pippippimpimpmiss_flatTree_bggen_2017_01-ver03.root"
  treeName     = "pippippimpimpmiss"
  inputData    = ROOT.RDataFrame(treeName, treeFileName).Define("TrackFound", TRACK_FOUND_FILTER)

  filterTopologies = {
    ""                                             : None,
    "__2#pi^{#plus}2#pi^{#minus}p"                 : '(ThrownTopology.GetString() == "2#pi^{#plus}2#pi^{#minus}p")',
    "__2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0}]" : '(ThrownTopology.GetString() == "2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0}]")',
    "__bkg"                                        : '(ThrownTopology.GetString() != "2#pi^{#plus}2#pi^{#minus}p")'
  }
  for suffix, filter in filterTopologies.items():
    overlayCases(inputData, "TruthDeltaP",      axisTitles = "#it{p}^{miss}_{truth} #minus #it{p}^{miss}_{kin. fit} (GeV/c)",                      binning = (600, -6, 6),     additionalFilter = filter, pdfFileNameSuffix = suffix)
    overlayCases(inputData, "TruthDeltaPOverP", axisTitles = "(#it{p}^{miss}_{truth} #minus #it{p}^{miss}_{kin. fit}) / #it{p}^{miss}_{kin. fit}", binning = (500, -2, 2),     additionalFilter = filter, pdfFileNameSuffix = suffix)
    overlayCases(inputData, "TruthDeltaTheta",  axisTitles = "#it{#theta}^{miss}_{truth} #minus #it{#theta}^{miss}_{kin. fit} (deg)",              binning = (200, -100, 100), additionalFilter = filter, pdfFileNameSuffix = suffix)
    overlayCases(inputData, "TruthDeltaPhi",    axisTitles = "#it{#phi}^{miss}_{truth} #minus #it{#phi}^{miss}_{kin. fit} (deg)",                  binning = (360, -180, 180), additionalFilter = filter, pdfFileNameSuffix = suffix)

  plotTopologyHist(inputData, normalize = False)
  plotTopologyHist(inputData, normalize = True)

  cutsArgs = [
    {},  # no extra cut
    {"additionalFilter" : "(NmbUnusedShowers == 0)", "pdfFileNameSuffix" : "_noUnusedShowers"},  # no unused showers and hence no unused energy in calorimeters
    # # the two cuts below are equivalent to the one above
    # {"additionalFilter" : "(EnergyUnusedShowers == 0)", "pdfFileNameSuffix" : "_noEnergyUnusedShowers"},
    # {"additionalFilter" : "(NmbUnusedShowers == 0) and (EnergyUnusedShowers == 0)", "pdfFileNameSuffix" : "_noShowers"}
  ]
  for cutArgs in cutsArgs:
    overlayTopologies(inputData, "NmbUnusedShowers",            axisTitles = "Number of Unused Showers",                       binning = (11, -0.5, 10.5), **cutArgs)
    overlayTopologies(inputData, "EnergyUnusedShowers",         axisTitles = "Unused Shower Energy (GeV)",                     binning = (60, 0, 6),       **cutArgs)
    # overlayTopologies(inputData, "MissingMassSquared",          axisTitles = "(#it{m}^{miss}_{kin. fit})^{2} (GeV/c^{2})^{2}", binning = (50, -0.5, 4.5),  **cutArgs)
    overlayTopologies(inputData, "MissingMassSquared_Measured", axisTitles = "(#it{m}^{miss}_{measured})^{2} (GeV/c^{2})^{2}", binning = (50, -0.5, 4.5),  **cutArgs)
    overlayTopologies(inputData, "MissingMassSquared_Measured", axisTitles = "(#it{m}^{miss}_{measured})^{2} (GeV/c^{2})^{2}", binning = (50, -0.5, 4.5),  **cutArgs)

  sideBandArgs = {"weightVariable" : ("AccidWeightFactorSb", "1 - AccidWeightFactor"), "pdfFileNameSuffix" : "_Sb"}
  sideBandYTitle = "Number of Combos (RF-Sideband)"
  # plot1D(inputData, ("MissingMass", "sqrt(MissingMassSquared)"), axisTitles = "#it{m}^{miss}_{kin. fit} (GeV/c^{2})",                   binning = (100, 0, 2))
  # plot1D(inputData, ("MissingMass", "sqrt(MissingMassSquared)"), axisTitles = "#it{m}^{miss}_{kin. fit} (GeV/c^{2});" + sideBandYTitle, binning = (100, 0, 2), **sideBandArgs)
  # plot1D(inputData, "MissingMassSquared",                        axisTitles = "(#it{m}^{miss}_{kin. fit})^{2} (GeV/c^{2})^{2}",                   binning = (225, -0.5, 4))
  # plot1D(inputData, "MissingMassSquared",                        axisTitles = "(#it{m}^{miss}_{kin. fit})^{2} (GeV/c^{2})^{2};" + sideBandYTitle, binning = (225, -0.5, 4), **sideBandArgs)
  plot1D(inputData, ("MissingMass_Measured", "sqrt(MissingMassSquared_Measured)"), axisTitles = "#it{m}^{miss}_{measured} (GeV/c^{2})",                   binning = (100, 0, 2))
  plot1D(inputData, ("MissingMass_Measured", "sqrt(MissingMassSquared_Measured)"), axisTitles = "#it{m}^{miss}_{measured} (GeV/c^{2});" + sideBandYTitle, binning = (100, 0, 2), **sideBandArgs)
  plot1D(inputData, "MissingMassSquared_Measured",                                 axisTitles = "(#it{m}^{miss}_{measured})^{2} (GeV/c^{2})^{2}",                   binning = (225, -0.5, 4))
  plot1D(inputData, "MissingMassSquared_Measured",                                 axisTitles = "(#it{m}^{miss}_{measured})^{2} (GeV/c^{2})^{2};" + sideBandYTitle, binning = (225, -0.5, 4), **sideBandArgs)
  plot1D(inputData, "AccidWeightFactor",                                           axisTitles = "RF Weight", binning = (1000, -2, 2), weightVariable = None)

  beamEnergyRange   = (3.0, 12.0)  # [GeV]
  nmbBeamEnergyBins = 9            # 1 GeV bin width
  energyBinWidth    = (beamEnergyRange[1] - beamEnergyRange[0]) / float(nmbBeamEnergyBins)
  histDef           = {"variable" : "MissingMassSquared_Measured", "axisTitles" : "(#it{m}^{miss}_{measured})^{2} (GeV/c^{2})^{2}", "binning" : (250, -0.5, 4.5)}
  for case, caseFilter in FILTER_CASES.items():
    caseData = inputData.Filter(caseFilter)
    plot1D(caseData, **histDef, pdfFileNameSuffix = f"_{case}")
    # plot distributions for beam-energy bins
    for energyBin in range(nmbBeamEnergyBins):
      energyBinMin = beamEnergyRange[0] + energyBin * energyBinWidth
      energyBinMax = energyBinMin + energyBinWidth
      energyFilter = f"(({energyBinMin} < BeamEnergy) and (BeamEnergy < {energyBinMax}))"
      energyBinData = caseData.Filter(energyFilter)
      plot1D(energyBinData, **histDef, pdfFileNameSuffix = f"_Egamma_{energyBinMin}_{energyBinMax}_{case}")

  plot2D(inputData, xVariable = "MissingProtonTheta",          yVariable = "MissingProtonP",            axisTitles = "#it{#theta}^{miss}_{kin. fit} (deg);#it{p}^{miss}_{kin. fit} (GeV/c)",  binning = (180, 0, 90, 400, 0, 9))
  plot2D(inputData, xVariable = "MissingProtonTheta",          yVariable = "MissingProtonPhi",          axisTitles = "#it{#theta}^{miss}_{kin. fit} (deg);#it{#phi}^{miss}_{kin. fit} (deg)", binning = (180, 0, 90, 360, -180, 180))
  plot2D(inputData, xVariable = "MissingProtonTheta_Measured", yVariable = "MissingProtonP_Measured",   axisTitles = "#it{#theta}^{miss}_{measured} (deg);#it{p}^{miss}_{measured} (GeV/c)",  binning = (180, 0, 90, 400, 0, 9))
  plot2D(inputData, xVariable = "MissingProtonTheta_Measured", yVariable = "MissingProtonPhi_Measured", axisTitles = "#it{#theta}^{miss}_{measured} (deg);#it{#phi}^{miss}_{measured} (deg)", binning = (180, 0, 90, 360, -180, 180))

  plot1D(inputData, "UnusedDeltaP",      axisTitles = "#it{p}^{miss}_{unused} #minus #it{p}^{miss}_{kin. fit} (GeV/c)",                      binning = (600, -6, 6))
  plot1D(inputData, "UnusedDeltaPOverP", axisTitles = "(#it{p}^{miss}_{unused} #minus #it{p}^{miss}_{kin. fit}) / #it{p}^{miss}_{kin. fit}", binning = (500, -2, 2))
  plot1D(inputData, "UnusedDeltaTheta",  axisTitles = "#it{#theta}^{miss}_{unused} #minus #it{#theta}^{miss}_{kin. fit} (deg)",              binning = (200, -100, 100))
  plot1D(inputData, "UnusedDeltaPhi",    axisTitles = "#it{#phi}^{miss}_{unused} #minus #it{#phi}^{miss}_{kin. fit} (deg)",                  binning = (360, -180, 180))
  unusedTrackData = inputData.Filter("(NmbUnusedTracks == 1)")  # make sure unused track info exists; NOTE! this assumes that there is maximum 1 unused track
  plot2D(unusedTrackData, xVariable = ("UnusedP_",     "UnusedP[0]"),     yVariable = "MissingProtonP",     axisTitles = "#it{p}^{miss}_{unused} (GeV/c);#it{p}^{miss}_{kin. fit} (GeV/c)",       binning = (400, 0, 9, 400, 0, 9))
  plot2D(unusedTrackData, xVariable = ("UnusedTheta_", "UnusedTheta[0]"), yVariable = "MissingProtonTheta", axisTitles = "#it{#theta}^{miss}_{unused} (deg);#it{#theta}^{miss}_{kin. fit} (deg)", binning = (360, 0, 180, 360, 0, 180))
  plot2D(unusedTrackData, xVariable = ("UnusedPhi_",   "UnusedPhi[0]"),   yVariable = "MissingProtonPhi",   axisTitles = "#it{#phi}^{miss}_{unused} (deg);#it{#phi}^{miss}_{kin. fit} (deg)",     binning = (360, -180, 180, 360, -180, 180))

#TODO make 2D plots for measured theta and phi
