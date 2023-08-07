#!/usr/bin/env python3


from collections.abc import Iterable
import os
import subprocess
from typing import Any, Dict, Iterable, Tuple

import ROOT


# simplified versions of Paul's criteria
UNUSED_TRACK_FOUND_CONDITION: str = "(" \
  + "(NmbUnusedTracks == 1)" \
  + " and ((MissingProtonTheta < 5) or (abs(UnusedDeltaPhi[0]) <= 30))" \
  + " and (abs(UnusedDeltaTheta[0]) <= 30)" \
  + " and (abs(UnusedDeltaPOverP[0]) <= 0.6)" \
  + ")"
UNUSED_TRACK_FOUND_CONDITION_MEASURED: str = "(" \
  + "(NmbUnusedTracks == 1)" \
  + " and ((MissingProtonTheta_Measured < 5) or (abs(UnusedDeltaPhi_Measured[0]) <= 30))" \
  + " and (abs(UnusedDeltaTheta_Measured[0]) <= 30)" \
  + " and (abs(UnusedDeltaPOverP_Measured[0]) <= 0.6)" \
  + ")"

# filter expressions for track-found cases
FILTER_CASES: Dict[str, str] = {
  "Total"   : "(true)",
  "Found"   : "(TrackFound == true)",
  "Missing" : "(TrackFound == false)"
}

# black + 7 colorblind-friendly colors rom M. Okabe and K. Ito, "How to make figures and presentations that are friendly to color blind people," University of Tokyo, 2002.
# see also Bang Wong, https://www.nature.com/articles/nmeth.1618.pdf
#     https://davidmathlogic.com/colorblind
#     https://yoshke.org/blog/colorblind-friendly-diagrams
COLORS_CB_FRIENDLY: Tuple[str, ...] = (
  "#000000",  # black
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#009E73",  # bluish green
  "#CC79A7",  # reddish purple
  "#56B4E9",  # sky blue
  "#E69F00",  # orange
  "#F0E442",  # yellow
)

def getRootColor(hexColor: str) -> int:
  '''Returns ROOT color index for given hex string in form #RRGGBB; if color does not exist yet in ROOT it is created'''
  ROOT.TColor.SetColorThreshold(0)  # ensure GetColor() returns exact color  # type: ignore
  return ROOT.TColor.GetColor(hexColor)  # type: ignore

def getCbFriendlyRootColor(index: int) -> int:
  '''Returns ROOT color index for given index in colorblind-friendly palette'''
  return getRootColor(COLORS_CB_FRIENDLY[index])

# 11 filled marker styles; the float is a relative scaling factor to obtain equal apparent sizes
MARKERS_FILLED: Tuple[Tuple[int, float], ...] = (
  (ROOT.kFullCircle,            0.75),  # type: ignore
  (ROOT.kFullSquare,            0.70),  # type: ignore
  (ROOT.kFullDiamond,           1),     # type: ignore
  (ROOT.kFullCross,             0.85),  # type: ignore
  (ROOT.kFullCrossX,            0.85),  # type: ignore
  (ROOT.kFullStar,              1),     # type: ignore
  (ROOT.kFullFourTrianglesX,    0.90),  # type: ignore
  (ROOT.kFullFourTrianglesPlus, 0.85),  # type: ignore
  (ROOT.kFullTriangleUp,        0.85),  # type: ignore
  (ROOT.kFullTriangleDown,      0.85),  # type: ignore
  (ROOT.kFullDoubleDiamond,     1.10),  # type: ignore
)
# 11 open marker styles
MARKERS_OPEN: Tuple[Tuple[int, float], ...] = (
  (ROOT.kOpenCircle,            0.75),  # type: ignore
  (ROOT.kOpenSquare,            0.70),  # type: ignore
  (ROOT.kOpenDiamond,           1),     # type: ignore
  (ROOT.kOpenCross,             0.85),  # type: ignore
  (ROOT.kOpenCrossX,            0.85),  # type: ignore
  (ROOT.kOpenStar,              1),     # type: ignore
  (ROOT.kOpenFourTrianglesX,    0.90),  # type: ignore
  (ROOT.kOpenFourTrianglesPlus, 0.85),  # type: ignore
  (ROOT.kOpenTriangleUp,        0.85),  # type: ignore
  (ROOT.kOpenTriangleDown,      0.85),  # type: ignore
  (ROOT.kOpenDoubleDiamond,     1.10),  # type: ignore
)

def setCbFriendlyStyle(
  graphOrHist:   Any,
  styleIndex:    int,  # index that switches between styles
  skipBlack:     bool  = True,  # if set black color is not used
  markerSize:    float = 1.5,
  filledMarkers: bool  = True,
) -> None:
  '''Sets line color and marker style, color, and size of a TGraph or TH1 according to a style index'''
  nmbStyles = min(len(COLORS_CB_FRIENDLY) - (1 if skipBlack else 0), len(MARKERS_FILLED), len(MARKERS_OPEN))
  assert styleIndex < nmbStyles, f"The style index {styleIndex} goes beyond the maximum of {nmbStyles} styles that are implemented"
  graphOrHist.SetMarkerStyle(MARKERS_FILLED[styleIndex][0] if filledMarkers else MARKERS_OPEN[styleIndex][0])
  graphOrHist.SetMarkerSize(markerSize * MARKERS_FILLED[styleIndex][1] if filledMarkers else MARKERS_OPEN[styleIndex][1])
  color = getCbFriendlyRootColor(styleIndex + (1 if skipBlack else 0))
  graphOrHist.SetMarkerColor(color)
  graphOrHist.SetLineColor(color)


def printGitInfo() -> None:
  '''Prints directory of this file and git hash in this directory'''
  repoDir = os.path.dirname(os.path.abspath(__file__))
  gitInfo = subprocess.check_output(["git", "describe", "--always"], cwd = repoDir).strip().decode()
  print(f"Running code in '{repoDir}', git version '{gitInfo}'")


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
  inFileNames = ("pippippimpimpmiss.RD_2017_01-ver04_030730.root", "pippippimpimpmiss.MCbggen_2017_01-ver03.root")
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
    canv = ROOT.TCanvas("justin_Proton_4pi_mm2_MCbggen_overlay" + ("_" + case if case != "" else ""))
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
  weightVariable = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
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
    hist = getHistND(data, (variable,), setDefaultYAxisTitle(axisTitles), binning, weightVariable, FILTER_CASES[case],
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


# C++ helper functor to work around fact that RDataFrame cannot fill TObjSTring into histogram
# see https://sft.its.cern.ch/jira/browse/ROOT-10246
#TODO make this code thread-safe
# see https://root-forum.cern.ch/t/filling-histograms-in-parallel/35460/3
# and https://root.cern/doc/master/mt201__parallelHistoFill_8C.html
CPP_CODE = """
struct fillHistWithTObjString {

  fillHistWithTObjString(TH1& hist)
    : _hist(hist)
  { }

  void
  operator ()(const TObjString& s)
  {
    _hist.Fill(s.GetString().Data(), 1);
    return;
  }

  TH1& _hist;

};
"""
ROOT.gInterpreter.Declare(CPP_CODE)
CPP_CODE = """
struct fillHistWithTObjStringWeighted {

  fillHistWithTObjStringWeighted(TH1& hist)
    : _hist(hist)
  { }

  void
  operator ()(
    const TObjString& s,
    const Double_t    w
  ) {
    _hist.Fill(s.GetString().Data(), w);
    return;
  }

  TH1& _hist;

};
"""
ROOT.gInterpreter.Declare(CPP_CODE)


def getTopologyHist(
  inputData,  # RDataFrame
  weightVariable   = None,  # may be None (= no weighting), string with column name, or tuple with new column definition
  filterExpression = None,
  histNameSuffix   = ""
):
  #TODO the following does not work without recreating the RDataFrame
  # switch to single-threaded mode because the above code is not thread-safe (yet)
  nmbThreads = None
  if ROOT.IsImplicitMTEnabled():
    nmbThreads = ROOT.GetThreadPoolSize()
    ROOT.DisableImplicitMT()
  # apply additional filters, if defined
  data = inputData.Filter(filterExpression) if filterExpression else inputData
  if not isinstance(weightVariable, str) and isinstance(weightVariable, Iterable):
    # create new weight column
    data = data.Define(weightVariable[0], weightVariable[1])
  # create histogram
  variable = "ThrownTopology"
  histName = variable
  if isinstance(weightVariable, str):
    histName += f"_{weightVariable}"
  elif isinstance(weightVariable, Iterable):
    histName += f"_{weightVariable[0]}"
  if filterExpression:
    histName += f"_{filterExpression}"
  if histNameSuffix:
    histName += f"_{histNameSuffix}"
  hist = ROOT.TH1F(histName, "", 1, 0, 1)
  fillHistWithTObjString         = ROOT.fillHistWithTObjString        (hist)
  fillHistWithTObjStringWeighted = ROOT.fillHistWithTObjStringWeighted(hist)
  # fill histogram
  if not weightVariable:
    data.Foreach(fillHistWithTObjString, [variable])
  elif isinstance(weightVariable, str):
    # use existing weight column
    data.Foreach(fillHistWithTObjStringWeighted, [variable, weightVariable])
  elif isinstance(weightVariable, Iterable):
    # use new weight column
    data.Foreach(fillHistWithTObjStringWeighted, [variable, weightVariable[0]])
  hist.LabelsDeflate("X")
  hist.LabelsOption(">", "X")  # sort topologies by number od combos
  # get ordered list of topology names
  xAxis = hist.GetXaxis()
  topoNames = [xAxis.GetBinLabel(binIndex) for binIndex in range(1, xAxis.GetNbins() + 1)]
  # restore multithreading if it was enabled
  if nmbThreads:
    ROOT.EnableImplicitMT(nmbThreads)
  return (topoNames, hist)


# helper function that returns dict { bin label : bin content }
def getCategorialTH1AsDict(hist):
  xAxis = hist.GetXaxis()
  return {xAxis.GetBinLabel(binIndex) : hist.GetBinContent(binIndex)
    for binIndex in range(1, xAxis.GetNbins() + 1)}


def plotTopologyHist(
  inputData,  # RDataFrame
  normalize         = False,
  maxNmbTopologies  = 10,
  additionalFilter  = None,
  pdfFileNamePrefix = "justin_Proton_4pi_",
  pdfFileNameSuffix = ""
):
  colorCases = {
    "Total"   : ROOT.kGray,
    "Found"   : ROOT.kGreen + 2,
    "Missing" : ROOT.kRed + 1
  }
  # get histogram data
  topoNames = {}  # dictionary of ordered list of topology names { case : [ topologyName ] }
  topoHists = {}  # dictionary of topology histograms { case : topologyHist }
  for case in FILTER_CASES.keys():
    caseData = inputData.Filter(FILTER_CASES[case])
    topoNames[case], topoHists[case] = getTopologyHist(caseData, weightVariable = "AccidWeightFactor", filterExpression = additionalFilter, histNameSuffix = case + ("_norm" if normalize else ""))
  # overlay distributions for cases
  hStack = ROOT.THStack(f"topologies",  ";;" + ("Fraction" if normalize else "Number") + " of Combos (RF-subtracted)" + (" [%]" if normalize else ""))
  topoLabels = topoNames["Total"]
  hists = {}  # memorize plots to print
  for case in FILTER_CASES.keys():
    # ensure that bin labels in all histograms have same order as defined by the "Total" histogram
    hist = ROOT.TH1F(f"{pdfFileNamePrefix}topologies_{case}{'_norm' if normalize else ''}{pdfFileNameSuffix}", case, len(topoLabels), 0, len(topoLabels))
    xAxis = hist.GetXaxis()
    for binIndex, binLabel in enumerate(topoLabels):
      xAxis.SetBinLabel(binIndex + 1, binLabel)
    # get histogram values as dictionary, i.e. { topology : count }
    histValues = getCategorialTH1AsDict(topoHists[case])
    # set bin content of histogram
    for binLabel in topoLabels:
      hist.Fill(binLabel, histValues[binLabel] if binLabel in histValues else 0)
    if normalize:
      hist.Scale(100 / hist.Integral())
    hist.SetLineColor(colorCases[case])
    if case == "Total":
      hist.SetFillColor(colorCases[case])
    hists[case] = hist
    hStack.Add(hist)
    print(f"plotTopologyHist(): {case} signal: {hist.GetBinContent(1)}{'%' if normalize else ' combos'}")
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}topologies{'_norm' if normalize else ''}{pdfFileNameSuffix}")
  hStack.Draw("NOSTACK HIST")
  hStack.SetMinimum(0)
  hStack.GetXaxis().SetRangeUser(0, maxNmbTopologies)
  # add legend
  legend = canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
  # add labels that show number or fraction outside of plot range
  legend.AddEntry(ROOT.MakeNullPointer(ROOT.TObject), "Other topologies:", "")
  for case in FILTER_CASES.keys():
    integralOtherTopos = hists[case].Integral(maxNmbTopologies, hists[case].GetNbinsX())
    legendEntry = legend.AddEntry(ROOT.MakeNullPointer(ROOT.TObject), "    " + str(round(integralOtherTopos)) + ("%" if normalize else " Combos"), "")
    legendEntry.SetTextColor(ROOT.kBlack if case == "Total" else colorCases[case])
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
  pdfFileNameSuffix = "_MCbggen_topologies"
):
  data = inputData.Filter(additionalFilter) if additionalFilter else inputData
  for case in FILTER_CASES.keys():
    caseData = data.Filter(FILTER_CASES[case])
    # get topologies with largest number of combos for given case
    #TODO fix call
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
  #TODO add command-line interface
  ROOT.gROOT.SetBatch(True)
  #TODO cannot change multithreading for existing data frame
  # ROOT.EnableImplicitMT(20)  # activate implicit multi-threading for RDataFrame; disable using ROOT.DisableImplicitMT()
  setupPlotStyle()

  # overlayMissingMassSquared()

  # dataset = None
  # dataset = "MCbggen_2017_01-ver03"
  # dataset = "MCbggen_2018_01-ver02"
  # isMonteCarlo = isMcBggen = True
  # dataset = "RD_2017_01-ver04_030730"
  # dataset = "RD_2018_01-ver02_041003"
  dataset = "RD_2019_11-ver01_071592"
  isMonteCarlo = isMcBggen = False
  histFileName = f"pippippimpimpmiss.{dataset}.root"          if dataset else "pippippimpimpmiss.root"
  treeFileName = f"pippippimpimpmiss_flatTree.{dataset}.root" if dataset else "pippippimpimpmiss_flatTree.root"
  treeName     = "pippippimpimpmiss"
  inputData    = ROOT.RDataFrame(treeName, treeFileName) \
                     .Define("TrackFound", UNUSED_TRACK_FOUND_CONDITION) \
                     .Filter("(-0.25 < MissingMassSquared_Measured) and (MissingMassSquared_Measured < 3.75)")  # limit data to fit range

  #TODO determine isMonteCarlo and isMcBggen flags from data
  if isMonteCarlo:
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

    if isMcBggen:

      cutsArgs = [
        {},  # no extra cut
        {"additionalFilter" : "(NmbUnusedShowers == 0)", "pdfFileNameSuffix" : "_noUnusedShowers"},  # no unused showers and hence no unused energy in calorimeters
      ]
      for cutArgs in cutsArgs:
        plotTopologyHist(inputData, normalize = False, **cutArgs)
        plotTopologyHist(inputData, normalize = True,  **cutArgs)

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
  plot1D(inputData, "AccidWeightFactor",                                           axisTitles = "RF Weight", binning = (1000, -2, 2), weightVariable = None)
  # plot1D(inputData, ("MissingMass", "sqrt(MissingMassSquared)"), axisTitles = "#it{m}^{miss}_{kin. fit} (GeV/c^{2})",                   binning = (100, 0, 2))
  # plot1D(inputData, ("MissingMass", "sqrt(MissingMassSquared)"), axisTitles = "#it{m}^{miss}_{kin. fit} (GeV/c^{2});" + sideBandYTitle, binning = (100, 0, 2), **sideBandArgs)
  # plot1D(inputData, "MissingMassSquared",                        axisTitles = "(#it{m}^{miss}_{kin. fit})^{2} (GeV/c^{2})^{2}",                   binning = (225, -0.5, 4))
  # plot1D(inputData, "MissingMassSquared",                        axisTitles = "(#it{m}^{miss}_{kin. fit})^{2} (GeV/c^{2})^{2};" + sideBandYTitle, binning = (225, -0.5, 4), **sideBandArgs)
  plot1D(inputData, ("MissingMass_Measured", "sqrt(MissingMassSquared_Measured)"), axisTitles = "#it{m}^{miss}_{measured} (GeV/c^{2})",                   binning = (100, 0, 2))
  plot1D(inputData, ("MissingMass_Measured", "sqrt(MissingMassSquared_Measured)"), axisTitles = "#it{m}^{miss}_{measured} (GeV/c^{2});" + sideBandYTitle, binning = (100, 0, 2), **sideBandArgs)

  # missing-mass squared distributions
  mm2HistDef         = {"variable" : "MissingMassSquared_Measured", "axisTitles" : "(#it{m}^{miss}_{measured})^{2} (GeV/c^{2})^{2}",                   "binning" : (125, -0.5, 4.5)}
  mm2HistDefSideBand = {"variable" : "MissingMassSquared_Measured", "axisTitles" : "(#it{m}^{miss}_{measured})^{2} (GeV/c^{2})^{2};" + sideBandYTitle, "binning" : (125, -0.5, 4.5), "weightVariable" : ("AccidWeightFactorSb", "1 - AccidWeightFactor")}
  overlayCases(inputData, **mm2HistDef)
  overlayCases(inputData, **mm2HistDefSideBand, pdfFileNameSuffix = f"_Sb")
  # plot overall distributions for each case
  for case, caseFilter in FILTER_CASES.items():
    caseData = inputData.Filter(caseFilter)
    plot1D(caseData, **mm2HistDef,         pdfFileNameSuffix = f"_{case}")
    plot1D(caseData, **mm2HistDefSideBand, pdfFileNameSuffix = f"_{case}_Sb")
  kinematicBinnings  = [
    # beam energy
    # {"variable" : "BeamEnergy",         "label" : "Beam Energy",                   "unit" : "GeV",   "nmbBins" :  9, "range" : (3.0, 12.0)},  # spring 2017
    {"variable" : "BeamEnergy",         "label" : "Beam Energy",                   "unit" : "GeV",   "nmbBins" : 10, "range" : (5.5, 11.5)},  # spring 2018
    # momentum of missing proton
    {"variable" : "MissingProtonP",     "label" : "#it{p}^{miss}_{kin. fit}",      "unit" : "GeV/c", "nmbBins" : 10, "range" : (0, 3.5)},
    # polar angle of missing proton
    {"variable" : "MissingProtonTheta", "label" : "#it{#theta}^{miss}_{kin. fit}", "unit" : "deg",   "nmbBins" : 10, "range" : (0, 65)},
    # azimuthal angle of missing proton
    {"variable" : "MissingProtonPhi",   "label" : "#it{#phi}^{miss}_{kin. fit}",   "unit" : "deg",   "nmbBins" : 10, "range" : (-180, +180)},
  ]
  for kinematicBinning in kinematicBinnings:
    kinBinVariable = kinematicBinning["variable"]
    nmbKinBins     = kinematicBinning["nmbBins"]
    kinBinRange    = kinematicBinning["range"]
    kinBinWidth    = (kinBinRange[1] - kinBinRange[0]) / float(nmbKinBins)
    # plot distributions for kinematic bins
    for kinBinIndex in range(nmbKinBins):
      kinBinMin = kinBinRange[0] + kinBinIndex * kinBinWidth
      kinBinMax = kinBinMin + kinBinWidth
      kinBinFilter = f"(({kinBinMin} < {kinBinVariable}) and ({kinBinVariable} < {kinBinMax}))"
      kinBinData = caseData.Filter(kinBinFilter)
      overlayCases(kinBinData, **mm2HistDef, pdfFileNameSuffix = f"_{kinBinVariable}_{kinBinMin}_{kinBinMax}")

  plot2D(inputData, xVariable = "MissingProtonTheta",          yVariable = "MissingProtonP",            axisTitles = "#it{#theta}^{miss}_{kin. fit} (deg);#it{p}^{miss}_{kin. fit} (GeV/c)",  binning = (180, 0, 90, 400, 0, 9))
  plot2D(inputData, xVariable = "MissingProtonTheta",          yVariable = "MissingProtonPhi",          axisTitles = "#it{#theta}^{miss}_{kin. fit} (deg);#it{#phi}^{miss}_{kin. fit} (deg)", binning = (180, 0, 90, 360, -180, 180))
  plot2D(inputData, xVariable = "MissingProtonTheta_Measured", yVariable = "MissingProtonP_Measured",   axisTitles = "#it{#theta}^{miss}_{measured} (deg);#it{p}^{miss}_{measured} (GeV/c)",  binning = (180, 0, 90, 400, 0, 9))
  plot2D(inputData, xVariable = "MissingProtonTheta_Measured", yVariable = "MissingProtonPhi_Measured", axisTitles = "#it{#theta}^{miss}_{measured} (deg);#it{#phi}^{miss}_{measured} (deg)", binning = (180, 0, 90, 360, -180, 180))

  plot1D(inputData, "MissingProtonP",      axisTitles = "#it{p}^{miss}_{kin. fit} (GeV/c)",    binning = (1000, 0, 10),  additionalFilter = "(NmbUnusedShowers == 0)")
  plot1D(inputData, "MissingProtonTheta",  axisTitles = "#it{#theta}^{miss}_{kin. fit} (deg)", binning = (1000, 0, 100), additionalFilter = "(NmbUnusedShowers == 0)")

  plot1D(inputData, "UnusedDeltaP",      axisTitles = "#it{p}^{miss}_{unused} #minus #it{p}^{miss}_{kin. fit} (GeV/c)",                      binning = (600, -6, 6))
  plot1D(inputData, "UnusedDeltaPOverP", axisTitles = "(#it{p}^{miss}_{unused} #minus #it{p}^{miss}_{kin. fit}) / #it{p}^{miss}_{kin. fit}", binning = (500, -2, 2))
  plot1D(inputData, "UnusedDeltaTheta",  axisTitles = "#it{#theta}^{miss}_{unused} #minus #it{#theta}^{miss}_{kin. fit} (deg)",              binning = (200, -100, 100))
  plot1D(inputData, "UnusedDeltaPhi",    axisTitles = "#it{#phi}^{miss}_{unused} #minus #it{#phi}^{miss}_{kin. fit} (deg)",                  binning = (360, -180, 180))
  overlayCases(inputData, "UnusedDeltaP",      axisTitles = "#it{p}^{miss}_{unused} #minus #it{p}^{miss}_{kin. fit} (GeV/c)",                      binning = (600, -6, 6))
  overlayCases(inputData, "UnusedDeltaPOverP", axisTitles = "(#it{p}^{miss}_{unused} #minus #it{p}^{miss}_{kin. fit}) / #it{p}^{miss}_{kin. fit}", binning = (500, -2, 2))
  overlayCases(inputData, "UnusedDeltaTheta",  axisTitles = "#it{#theta}^{miss}_{unused} #minus #it{#theta}^{miss}_{kin. fit} (deg)",              binning = (200, -100, 100))
  overlayCases(inputData, "UnusedDeltaPhi",    axisTitles = "#it{#phi}^{miss}_{unused} #minus #it{#phi}^{miss}_{kin. fit} (deg)",                  binning = (360, -180, 180))
  unusedTrackData = inputData.Filter("(NmbUnusedTracks == 1)")  # make sure unused track info exists; NOTE! this assumes that there is maximum 1 unused track
  plot2D(unusedTrackData, xVariable = ("UnusedP_",     "UnusedP[0]"),     yVariable = "MissingProtonP",     axisTitles = "#it{p}^{miss}_{unused} (GeV/c);#it{p}^{miss}_{kin. fit} (GeV/c)",       binning = (400, 0, 9, 400, 0, 9))
  plot2D(unusedTrackData, xVariable = ("UnusedTheta_", "UnusedTheta[0]"), yVariable = "MissingProtonTheta", axisTitles = "#it{#theta}^{miss}_{unused} (deg);#it{#theta}^{miss}_{kin. fit} (deg)", binning = (360, 0, 180, 360, 0, 180))
  plot2D(unusedTrackData, xVariable = ("UnusedPhi_",   "UnusedPhi[0]"),   yVariable = "MissingProtonPhi",   axisTitles = "#it{#phi}^{miss}_{unused} (deg);#it{#phi}^{miss}_{kin. fit} (deg)",     binning = (360, -180, 180, 360, -180, 180))

#TODO make 2D plots for measured theta and phi
