#!/usr/bin/env python3


from collections.abc import Sequence
import ctypes
import functools
import os
import subprocess
from typing import (
  Any,
  Dict,
  List,
  Optional,
  Sequence,
  Tuple,
  Union,
)

import ROOT


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


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
  "Missing" : "(TrackFound == false)",
}

COLOR_CASES = {
  "Total"   : ROOT.kGray,
  "Found"   : ROOT.kGreen + 2,
  "Missing" : ROOT.kRed + 1,
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
  """Returns ROOT color index for given hex string in form #RRGGBB; if color does not exist yet in ROOT it is created"""
  ROOT.TColor.SetColorThreshold(0)  # ensure GetColor() returns exact color
  return ROOT.TColor.GetColor(hexColor)

def getCbFriendlyRootColor(index: int) -> int:
  """Returns ROOT color index for given index in colorblind-friendly palette"""
  return getRootColor(COLORS_CB_FRIENDLY[index])

# 11 filled marker styles; the float is a relative scaling factor to obtain equal apparent sizes
MARKERS_FILLED: Tuple[Tuple[int, float], ...] = (
  (ROOT.kFullCircle,            0.75),
  (ROOT.kFullSquare,            0.70),
  (ROOT.kFullDiamond,           1.00),
  (ROOT.kFullCross,             0.85),
  (ROOT.kFullCrossX,            0.85),
  (ROOT.kFullStar,              1.00),
  (ROOT.kFullFourTrianglesX,    0.90),
  (ROOT.kFullFourTrianglesPlus, 0.85),
  (ROOT.kFullTriangleUp,        0.85),
  (ROOT.kFullTriangleDown,      0.85),
  (ROOT.kFullDoubleDiamond,     1.10),
)
# 11 open marker styles
MARKERS_OPEN: Tuple[Tuple[int, float], ...] = (
  (ROOT.kOpenCircle,            0.75),
  (ROOT.kOpenSquare,            0.70),
  (ROOT.kOpenDiamond,           1.00),
  (ROOT.kOpenCross,             0.85),
  (ROOT.kOpenCrossX,            0.85),
  (ROOT.kOpenStar,              1.00),
  (ROOT.kOpenFourTrianglesX,    0.90),
  (ROOT.kOpenFourTrianglesPlus, 0.85),
  (ROOT.kOpenTriangleUp,        0.85),
  (ROOT.kOpenTriangleDown,      0.85),
  (ROOT.kOpenDoubleDiamond,     1.10),
)

def setCbFriendlyStyle(
  graphOrHist:   Any,
  styleIndex:    int,  # index that switches between styles
  skipBlack:     bool  = True,  # if set black color is not used
  markerSize:    float = 1.5,
  filledMarkers: bool  = True,
) -> None:
  """Sets line color and marker style, color, and size of a TGraph or TH1 according to a style index"""
  nmbStyles = min(len(COLORS_CB_FRIENDLY) - (1 if skipBlack else 0), len(MARKERS_FILLED), len(MARKERS_OPEN))
  assert styleIndex < nmbStyles, f"The style index {styleIndex} goes beyond the maximum of {nmbStyles} styles that are implemented"
  graphOrHist.SetMarkerStyle(MARKERS_FILLED[styleIndex][0] if filledMarkers else MARKERS_OPEN[styleIndex][0])
  graphOrHist.SetMarkerSize(markerSize * MARKERS_FILLED[styleIndex][1] if filledMarkers else MARKERS_OPEN[styleIndex][1])
  color = getCbFriendlyRootColor(styleIndex + (1 if skipBlack else 0))
  graphOrHist.SetMarkerColor(color)
  graphOrHist.SetLineColor(color)


def printGitInfo() -> None:
  """Prints directory of this file and git hash in this directory"""
  repoDir = os.path.dirname(os.path.abspath(__file__))
  gitInfo = subprocess.check_output(["git", "describe", "--always"], cwd = repoDir).strip().decode()
  print(f"Running code in '{repoDir}', git version '{gitInfo}'")


def setupPlotStyle() -> None:
  """Defines ROOT plot style"""
  #TODO remove dependency from external file or add file to repo
  ROOT.gROOT.LoadMacro("./rootlogon.C")
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


def drawZeroLine(obj, style = ROOT.kDashed, color = ROOT.kBlack) -> None:
  """Draws zero line when necessary"""
  objType = obj.IsA().GetName()
  if (objType == "TCanvas") or (objType == "TPad"):
    xMin = ctypes.c_double()
    xMax = ctypes.c_double()
    yMin = ctypes.c_double()
    yMax = ctypes.c_double()
    obj.GetRangeAxis(xMin, yMin, xMax, yMax)
    if (yMin.value < 0) and (yMax.value > 0):
      zeroLine = ROOT.TLine()
      zeroLine.SetLineStyle(style)
      zeroLine.SetLineColor(color)
      return zeroLine.DrawLine(xMin, 0, xMax, 0)
  elif objType.startswith("TH") or objType.startswith("TGraph") or objType.startswith("TMulti"):
    xAxis = obj.GetXaxis()
    yAxis = obj.GetYaxis()
    if (yAxis.GetXmin() < 0) and (yAxis.GetXmax() > 0):
      zeroLine = ROOT.TLine()
      zeroLine.SetLineStyle(style)
      zeroLine.SetLineColor(color)
      return zeroLine.DrawLine(xAxis.GetBinLowEdge(xAxis.GetFirst()), 0, xAxis.GetBinUpEdge(xAxis.GetLast()), 0)
    elif (yAxis.GetXmin() > 0) and (yAxis.GetXmax() > 0):
      obj.SetMinimum(0)
  else:
    raise TypeError(f"drawZeroLine() not (yet) implemented for object of type '{objType}'")


def callMemberFunctionsWithArgs(
  instance:          Any,             # instance for which member functions will be called
  functionsWithArgs: Dict[str, Any],  # member-function names with argments
) -> None:
  """Calls member functions of given object with given arguments"""
  for functionName, argument in functionsWithArgs.items():
    function = getattr(instance, functionName, None)
    if function is None:
      continue
    # print(f"Calling member function '{functionName}({argument})' of {instance}")
    function(argument)


def overlayDataSamples1D(
  dataSamples:       Dict[str, Dict[str, Any]],  # file name and style definition for each data-set label
  treeName:          str,  # tree name to read
  variable:          Union[str, Tuple[str, str]],  # variable to plot; may be column name, or tuple with new column definition
  axisTitles:        str,  # semicolon-separated list
  binning:           Tuple[int, float, float],  # tuple with binning definition
  weightVariable:    Optional[Union[str, Tuple[str, str]]] = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
  pdfFileNamePrefix: str = "justin_Proton_4pi_",
  pdfFileNameSuffix: str = "",
  additionalFilter:  Optional[str] = None,
) -> None:
  """Overlays 1D histograms from given data samples"""
  print(f"Overlaying distributions for '{variable}' for data samples {', '.join(dataSamples.keys())}")
  hStack = ROOT.THStack(f"{variable}", ";" + setDefaultYAxisTitle(axisTitles))
  hists: List[ROOT.TH1D] = []  # keep histograms in memory
  normIntegral = None  # index of histogram to normalize to
  for dataLabel, dataSample in dataSamples.items():
    # get tree
    treeFileName = dataSample["fileName"]
    print(f"Reading data sample '{dataLabel}' from tree '{treeName}' in file '{treeFileName}'")
    data = ROOT.RDataFrame(treeName, treeFileName) \
               .Define("TrackFound", UNUSED_TRACK_FOUND_CONDITION) \
               .Filter("(-0.25 < MissingMassSquared_Measured) and (MissingMassSquared_Measured < 3.75)")  # limit data to fit range
    hist: ROOT.TH1D = getHistND(data, (variable,), setDefaultYAxisTitle(axisTitles), binning, weightVariable, additionalFilter,
                                histNameSuffix = dataLabel, histTitle = dataLabel)
    callMemberFunctionsWithArgs(hist, dataSample)
    if dataSample.get("normToThis", None) is not None:
      print(f"Normalizing all histograms to '{dataLabel}'")
      normIntegral = hist.Integral()
    hists.append(hist)
    hStack.Add(hist.GetPtr())
  # normalize histograms
  if normIntegral is not None:
    for hist in hists:
      hist.Scale(normIntegral / hist.Integral())
  # draw distributions
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{variable}_overlay_{'_'.join(dataSamples.keys())}{pdfFileNameSuffix}")
  hStack.Draw("NOSTACK HIST")
  # add legend
  # canv.BuildLegend()  # automatic placement with width 0.3 and height 0.21
  canv.BuildLegend(0.3, 0.15, 0.3, 0.15)  # automatic placement with width 0.3 and height 0.15
  # canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
  canv.SaveAs(".pdf")


def drawHistogram(
  inFileName:        str,
  histName:          str,
  rebinFactor:       Union[int, Sequence[int]] = 1,  # if integer -> rebin x axis; if sequence of 2 integers -> rebin x and y axes
  drawOption:        str = "HIST",
  pdfFileNamePrefix: str = "justin_Proton_4pi_",
  pdfFileNameSuffix: str = "",
) -> None:
  """Plots histogram with given name in ROOT file with given name"""
  # get histogram
  inFile = ROOT.TFile(inFileName)
  hist = inFile.Get(histName)
  if isinstance(hist, ROOT.TH2):
    if isinstance(rebinFactor, int):
      hist.RebinX(rebinFactor)
    elif isinstance(rebinFactor, Sequence):
      hist.Rebin2D(rebinFactor[0], rebinFactor[1])
  elif isinstance(hist, ROOT.TH1):
    hist.Rebin(rebinFactor)
  # draw histogram
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{hist.GetName()}{pdfFileNameSuffix}")
  hist.Draw(drawOption)
  canv.SaveAs(".pdf")


def getHistND(
  inputData:        ROOT.RDataFrame,
  variables:        Tuple[Union[str, Tuple[str, str]], ...],  # variable(s) to plot; is tuple of either column names or tuples with new column definitions; defines dimension of histogram
  axisTitles:       str,    # semicolon-separated list
  binning:          Tuple,  # tuple with binning definitions
  weightVariable:   Optional[Union[str, Tuple[str, str]]] = None,  # may be None (= no weighting), string with column name, or tuple with new column definition
  filterExpression: Optional[str] = None,
  histNameSuffix:   str = "",
  histTitle:        str = "",
) -> Union[ROOT.TH1D, ROOT.TH2D]:
  """Creates histogram from given variables in RDataFrame, applying optional weighting and filtering"""
  histDim = len(variables)
  assert 1 <= histDim <= 2, "currently, only 1D and 2D histograms are supported"
  # apply additional filters, if defined
  data = inputData.Filter(filterExpression) if filterExpression else inputData
  columnNames: List[str] = [""] * len(variables)
  for index, variable in enumerate(variables):
    if isinstance(variable, str):
      # use existing variable column
      columnNames[index] = variable
    elif isinstance(variable, Sequence):
      # create new variable column
      data  = data.Define(variable[0], variable[1])
      columnNames[index] = variable[0]
    assert columnNames[index] != "", f"failed to get column name for variable '{variable}'"
  if isinstance(weightVariable, Sequence) and not isinstance(weightVariable, str):
    # create new weight column
    data = data.Define(weightVariable[0], weightVariable[1])
  # create histogram
  hist = None
  histDef = (  # histogram definition
    "_".join(columnNames + ([histNameSuffix] if histNameSuffix else [])),
    f"{histTitle};{axisTitles}",
    *binning,
  )
  # get member function to create histogram
  HistoND = getattr(data, "Histo1D") if histDim == 1 else getattr(data, "Histo2D")
  if not weightVariable:
    hist = HistoND(histDef, *columnNames)
  elif isinstance(weightVariable, str):
    # use existing weight column
    hist = HistoND(histDef, *columnNames, weightVariable)
  elif isinstance(weightVariable, Sequence):
    # use new weight column
    hist = HistoND(histDef, *columnNames, weightVariable[0])
  assert hist is not None, f"failed to create histogram for weight variable '{weightVariable}'"
  return hist


def setDefaultYAxisTitle(
  axisTitles:    str,  # semicolon-separated list
  defaultYTitle: str = "Number of Combos (RF-subtracted)",
):
  """Sets default y-axis title if not provided by `axisTitles`"""
  titles = axisTitles.split(";")
  if (len(titles) == 1):
    return titles[0] + ";" + defaultYTitle
  elif (len(titles) == 0):
    return ";" + defaultYTitle
  else:
    return axisTitles


def plot1D(
  inputData:         ROOT.RDataFrame,
  variable:          Union[str, Tuple[str, str]],  # variable to plot; may be column name, or tuple with new column definition
  axisTitles:        str,  # semicolon-separated list
  binning:           Tuple[int, float, float],  # tuple with binning definition
  weightVariable:    Optional[Union[str, Tuple[str, str]]] = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
  pdfFileNamePrefix: str = "justin_Proton_4pi_",
  pdfFileNameSuffix: str = "",
  additionalFilter:  Optional[str] = None,
) -> None:
  """Plots 1D distribution for given variable, applying optional weighting and filtering"""
  hist: ROOT.TH1D = getHistND(inputData, (variable,), setDefaultYAxisTitle(axisTitles), binning, weightVariable, additionalFilter)
  # draw distributions
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{hist.GetName()}{pdfFileNameSuffix}")
  hist.Draw("HIST")
  drawZeroLine(hist)
  canv.SaveAs(".pdf")


def plot2D(
  inputData:         ROOT.RDataFrame,
  xVariable:         Union[str, Tuple[str, str]],  # x variable to plot; may be column name, or tuple with new column definition
  yVariable:         Union[str, Tuple[str, str]],  # y variable to plot; may be column name, or tuple with new column definition
  axisTitles:        str,  # semicolon-separated list
  binning:           Tuple[int, float, float, int, float, float],  # tuple with binning definition
  weightVariable:    Optional[Union[str, Tuple[str, str]]] = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
  pdfFileNamePrefix: str = "justin_Proton_4pi_",
  pdfFileNameSuffix: str = "",
  additionalFilter:  Optional[str] = None,
) -> None:
  """Plots 2D distribution for given x and y variables, applying optional weighting and filtering"""
  hist: ROOT.TH2D = getHistND(inputData, (xVariable, yVariable), axisTitles, binning, weightVariable, additionalFilter)
  # draw distributions
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{hist.GetName()}{pdfFileNameSuffix}")
  hist.Draw("COLZ")
  canv.SaveAs(".pdf")


def overlayCases(
  inputData:         ROOT.RDataFrame,
  variable:          Union[str, Tuple[str, str]],  # variable to plot; may be column name, or tuple with new column definition
  axisTitles:        str,  # semicolon-separated list
  binning:           Tuple[int, float, float],  # tuple with binning definition
  weightVariable:    Optional[Union[str, Tuple[str, str]]] = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
  pdfFileNamePrefix: str = "justin_Proton_4pi_",
  pdfFileNameSuffix: str = "",
  additionalFilter:  Optional[str] = None,
) -> None:
  """Overlays 1D distributions of given variable for "Total", "Found", and "Missing" cases"""
  data = inputData.Filter(additionalFilter) if additionalFilter else inputData
  # overlay distributions for cases
  hStack = ROOT.THStack(f"{variable}", ";" + setDefaultYAxisTitle(axisTitles))
  hists = []
  for case in FILTER_CASES.keys():
    hist: ROOT.TH1D = getHistND(data, (variable,), setDefaultYAxisTitle(axisTitles), binning, weightVariable, FILTER_CASES[case],
                                histNameSuffix = case, histTitle = case)
    hist.SetLineColor(COLOR_CASES[case])
    if case == "Total":
      hist.SetFillColor(COLOR_CASES[case])
    hists.append(hist)
    hStack.Add(hist.GetPtr())
  # draw distributions
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{variable}_cases{pdfFileNameSuffix}")
  hStack.Draw("NOSTACK HIST")
  # add legend
  canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
  drawZeroLine(hStack)
  canv.SaveAs(".pdf")


# C++ helper functors that fill TObjString into TH1
# workaround because RDataFrame cannot fill TObjSTring into histogram
# !Note! the code is not thread-safe; RDataFrame needs to run in single-threaded mode
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
  inputData:        ROOT.RDataFrame,
  weightVariable:   Optional[Union[str, Tuple[str, str]]] = None,  # may be None (= no weighting), string with column name, or tuple with new column definition
  filterExpression: Optional[str] = None,
  histNameSuffix:   str = "",
) -> Tuple[List[str], ROOT.TH1F]:
  """Fills categorical histogram with counts for each generated topology, applying optional weighting and filtering, and returns list of topology strings and histogram"""
  # apply additional filters, if defined
  data = inputData.Filter(filterExpression) if filterExpression else inputData
  if not isinstance(weightVariable, str) and isinstance(weightVariable, Sequence):
    # create new weight column
    data = data.Define(weightVariable[0], weightVariable[1])
  # create histogram
  variable = "ThrownTopology"
  histName = variable
  if isinstance(weightVariable, str):
    histName += f"_{weightVariable}"
  elif isinstance(weightVariable, Sequence):
    histName += f"_{weightVariable[0]}"
  if filterExpression:
    histName += f"_{filterExpression}"
  if histNameSuffix:
    histName += f"_{histNameSuffix}"
  hist = ROOT.TH1F(histName, "", 1, 0, 1)
  # fill histogram
  fillHistWithTObjString         = ROOT.fillHistWithTObjString        (hist)
  fillHistWithTObjStringWeighted = ROOT.fillHistWithTObjStringWeighted(hist)
  if not weightVariable:
    data.Foreach(fillHistWithTObjString, [variable])
  elif isinstance(weightVariable, str):
    # use existing weight column
    data.Foreach(fillHistWithTObjStringWeighted, [variable, weightVariable])
  elif isinstance(weightVariable, Sequence):
    # use new weight column
    data.Foreach(fillHistWithTObjStringWeighted, [variable, weightVariable[0]])
  hist.LabelsDeflate("X")
  hist.LabelsOption(">", "X")  # sort topologies by number od combos
  # get ordered list of topology names
  xAxis = hist.GetXaxis()
  topoNames = [xAxis.GetBinLabel(binIndex) for binIndex in range(1, xAxis.GetNbins() + 1)]
  return (topoNames, hist)


def getCategoricalTH1AsDict(hist: ROOT.TH1) -> Dict[str, float]:
  """Returns categorical histogram as dict { bin label : bin content }"""
  xAxis = hist.GetXaxis()
  return {xAxis.GetBinLabel(binIndex) : hist.GetBinContent(binIndex)
    for binIndex in range(1, xAxis.GetNbins() + 1)}


def plotTopologyHist(
  inputData:         ROOT.RDataFrame,
  normalize:         bool = False,
  maxNmbTopologies:  int  = 10,
  additionalFilter:  Optional[str] = None,
  pdfFileNamePrefix: str = "justin_Proton_4pi_",
  pdfFileNameSuffix: str = "",
) -> None:
  """Plots categorical histogram with counts or fraction for each generated topology, applying optional weighting and filtering"""
  # get histogram data
  topoNames: Dict[str, List[str]] = {}  # dictionary of ordered list of topology names { case : [ topologyName ] }
  topoHists: Dict[str, ROOT.TH1F] = {}  # dictionary of topology histograms { case : topologyHist }
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
    histValues = getCategoricalTH1AsDict(topoHists[case])
    # set bin content of histogram
    for binLabel in topoLabels:
      hist.Fill(binLabel, histValues[binLabel] if binLabel in histValues else 0)
    if normalize:
      hist.Scale(100 / hist.Integral())
    hist.SetLineColor(COLOR_CASES[case])
    if case == "Total":
      hist.SetFillColor(COLOR_CASES[case])
    hists[case] = hist
    hStack.Add(hist)
    print(f"plotTopologyHist(): {case} signal: {hist.GetBinContent(1)}{'%' if normalize else ' combos'}")
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}topologies{'_norm' if normalize else ''}{pdfFileNameSuffix}")
  hStack.Draw("NOSTACK HIST")
  hStack.SetMinimum(0)
  if normalize:
    hStack.SetMaximum(100)
  hStack.GetXaxis().SetRangeUser(0, maxNmbTopologies)
  # add legend
  legend = canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
  # add labels that show number or fraction outside of plot range
  legend.AddEntry(ROOT.MakeNullPointer(ROOT.TObject), "Other topologies:", "")
  for case in FILTER_CASES.keys():
    integralOtherTopos = hists[case].Integral(maxNmbTopologies, hists[case].GetNbinsX())
    legendEntry = legend.AddEntry(ROOT.MakeNullPointer(ROOT.TObject), "    " + str(round(integralOtherTopos)) + ("%" if normalize else " Combos"), "")
    legendEntry.SetTextColor(ROOT.kBlack if case == "Total" else COLOR_CASES[case])
  canv.SaveAs(".pdf")


def overlayTopologies(
  inputData:         ROOT.RDataFrame,
  variable:          Union[str, Tuple[str, str]],  # variable to plot; may be column name, or tuple with new column definition
  axisTitles:        str,  # semicolon-separated list
  binning:           Tuple[int, float, float],  # tuple with binning definition
  toposToPlot:       Dict[str, List[str]],  # topologies to plot for each case
  additionalFilter:  Optional[str] = None,
  pdfFileNamePrefix: str = "justin_Proton_4pi_",
  pdfFileNameSuffix: str = "_MCbggen_topologies",
) -> None:
  """Overlays 1D distributions for given variable from overall data sample and distributions for the `maxNmbTopologies` topologies with the largest number of combos from the bggen MC sample"""
  data = inputData.Filter(additionalFilter) if additionalFilter else inputData
  for case in FILTER_CASES.keys():
    caseData = data.Filter(FILTER_CASES[case])
    # get topologies with largest number of combos for given case
    hStack = ROOT.THStack(f"{variable}_{case}", f"{case};{setDefaultYAxisTitle(axisTitles)}")
    hists = []
    # overlay distributions for topologies
    for index, topo in enumerate(toposToPlot[case]):
      hist: ROOT.TH1D = getHistND(caseData, (variable,), setDefaultYAxisTitle(axisTitles), binning, "AccidWeightFactor",
                                  (f'ThrownTopology.GetString() == "{topo}"' if topo != "Total" else "true"), histNameSuffix = f"{case}_{topo}", histTitle = topo)
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
    drawZeroLine(hStack)
    canv.SaveAs(".pdf")


if __name__ == "__main__":
  #TODO add command-line interface
  ROOT.gROOT.SetBatch(True)
  #TODO cannot change multithreading for existing data frame
  # ROOT.EnableImplicitMT(20)  # activate implicit multi-threading for RDataFrame; disable using ROOT.DisableImplicitMT()
  setupPlotStyle()

  dataSamplesToOverlay = {
    "bggen MC (scaled)" : {
      "fileName" : "./data/MCbggen/2018_01-ver02/pippippimpimpmiss_flatTree.MCbggen_2018_01-ver02.root",
      # define plot style
      "SetLineColor" : ROOT.kGray,
      "SetFillColor" : ROOT.kGray,
    },
    "Real Data" : {
      "fileName"   : "./data/RD/2018_01-ver02/pippippimpimpmiss_flatTree.RD_2018_01-ver02.root",
      "normToThis" : True,
    },
  }
  overlayArgs = {
    "dataSamples" : dataSamplesToOverlay,
    "treeName"    : "pippippimpimpmiss",
  }
  overlayDataSamples1D(variable = "NmbUnusedShowers", axisTitles = "Number of Unused Showers", binning = (11, -0.5, 10.5), **overlayArgs)
  overlayArgs = {
    **overlayArgs,
    "additionalFilter"  : "(NmbUnusedShowers == 0)",
    "pdfFileNameSuffix" : "_noUnusedShowers",
  }
  overlayDataSamples1D(variable = "KinFitPVal",         axisTitles = "#it{#chi}^{2}_{kim. fit} #it{P}-value", binning = (150, 0, 1),      **overlayArgs)
  overlayDataSamples1D(variable = "MissingProtonP",     axisTitles = "#it{p}^{miss}_{kin. fit} (GeV/#it{c})", binning = (500, 0, 10),     **overlayArgs)
  overlayDataSamples1D(variable = "MissingProtonTheta", axisTitles = "#it{#theta}^{miss}_{kin. fit} (deg)",   binning = (200, 0, 100),    **overlayArgs)
  overlayDataSamples1D(variable = "MissingProtonPhi",   axisTitles = "#it{#phi}^{miss}_{kin. fit} (deg)",     binning = (180, -180, 180), **overlayArgs)
  for case, caseFilter in FILTER_CASES.items():
    overlayDataSamples1D(dataSamplesToOverlay, treeName = "pippippimpimpmiss", variable = "MissingMassSquared_Measured",
                         axisTitles = "(#it{m}^{miss}_{measured})^{2} (GeV/#it{c}^{2})^{2}", binning = (125, -0.5, 4.5),
                         additionalFilter = f"((NmbUnusedShowers == 0) and {caseFilter})", pdfFileNameSuffix = f"_{case}_noUnusedShowers")

  maxNmbTopologies = 10
  # dataSet = {}
  # dataset = {"type" : "MCbggen", "period" : "2017_01-ver03"}
  dataSet = {"type" : "MCbggen", "period" : "2018_01-ver02"}
  isMonteCarlo = isMcBggen = True
  # dataset = "RD_2017_01-ver04_030730"
  # dataset = "RD_2018_01-ver02_041003"
  # dataset = "RD_2019_11-ver01_071592"
  # dataSet = {"type" : "RD", "period" : "2018_01-ver02"}
  # isMonteCarlo = isMcBggen = False
  treeName     = "pippippimpimpmiss"
  treeFileName = f"./{treeName}_flatTree.root" if not dataSet else f"./data/{dataSet['type']}/{dataSet['period']}/{treeName}_flatTree.{dataSet['type']}_{dataSet['period']}.root"
  print(f"Reading tree '{treeName}' in file '{treeFileName}'")
  inputData    = ROOT.RDataFrame(treeName, treeFileName) \
                     .Define("TrackFound", UNUSED_TRACK_FOUND_CONDITION) \
                     .Filter("(-0.25 < MissingMassSquared_Measured) and (MissingMassSquared_Measured < 3.75)")  # limit data to fit range

  #TODO determine isMonteCarlo and isMcBggen flags from data
  if isMonteCarlo:
    filterTopologies = {
      ""                                             : None,
      # "__2#pi^{#plus}2#pi^{#minus}p"                 : '(ThrownTopology.GetString() == "2#pi^{#plus}2#pi^{#minus}p")',
      # "__2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0}]" : '(ThrownTopology.GetString() == "2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0}]")',
      "__bkg"                                        : '(ThrownTopology.GetString() != "2#pi^{#plus}2#pi^{#minus}p")',
    }
    for suffix, filter in filterTopologies.items():
      # overlayCases(inputData, "TruthDeltaP",      axisTitles = "#it{p}^{miss}_{truth} #minus #it{p}^{miss}_{kin. fit} (GeV/#it{c})",                 binning = (600, -6, 6),     additionalFilter = filter, pdfFileNameSuffix = suffix)
      overlayCases(inputData, "TruthDeltaPOverP", axisTitles = "(#it{p}^{miss}_{truth} #minus #it{p}^{miss}_{kin. fit}) / #it{p}^{miss}_{kin. fit}", binning = (500, -2, 2),     additionalFilter = filter, pdfFileNameSuffix = suffix)
      overlayCases(inputData, "TruthDeltaTheta",  axisTitles = "#it{#theta}^{miss}_{truth} #minus #it{#theta}^{miss}_{kin. fit} (deg)",              binning = (200, -100, 100), additionalFilter = filter, pdfFileNameSuffix = suffix)
      overlayCases(inputData, "TruthDeltaPhi",    axisTitles = "#it{#phi}^{miss}_{truth} #minus #it{#phi}^{miss}_{kin. fit} (deg)",                  binning = (360, -180, 180), additionalFilter = filter, pdfFileNameSuffix = suffix)

    if isMcBggen:

      cutsArgs: List[Dict[str, Any]] = [
        {},  # no extra cut
        {"additionalFilter" : "(NmbUnusedShowers == 0)", "pdfFileNameSuffix" : "_noUnusedShowers"},  # no unused showers and hence no unused energy in calorimeters
        {"additionalFilter" : "((NmbUnusedShowers == 0) and (BestMissingMatchDistTOF < 40))", "pdfFileNameSuffix" : "_noUnusedShowersMatchToF"},  # no unused showers and ToF hit within certain distance
      ]
      for kwargs in cutsArgs:
        plotTopologyHist(inputData, normalize = False, **kwargs)
        plotTopologyHist(inputData, normalize = True,  **kwargs)

      cutsArgs = [
        {},  # no extra cut
        {"additionalFilter" : "(NmbUnusedShowers == 0)", "pdfFileNameSuffix" : "_noUnusedShowers"},  # no unused showers and hence no unused energy in calorimeters
        {"additionalFilter" : "((NmbUnusedShowers == 0) and (BestMissingMatchDistTOF < 40))", "pdfFileNameSuffix" : "_noUnusedShowersMatchToF"},  # no unused showers and ToF hit within certain distance
        # # the two cuts below are equivalent to the one above
        # {"additionalFilter" : "(EnergyUnusedShowers == 0)", "pdfFileNameSuffix" : "_noEnergyUnusedShowers"},
        # {"additionalFilter" : "(NmbUnusedShowers == 0) and (EnergyUnusedShowers == 0)", "pdfFileNameSuffix" : "_noShowers"}
      ]
      for kwargs in cutsArgs:
        # get topologies with the largest number of combos for given case
        toposToPlot: Dict[str, List[str]] = {}
        for case in FILTER_CASES.keys():
          caseData = inputData.Filter(FILTER_CASES[case])
          toposToPlot[case], _ = getTopologyHist(caseData, filterExpression = kwargs.get("additionalFilter", None))
          toposToPlot[case] = ["Total"] + toposToPlot[case][:maxNmbTopologies]
        overlayTopologies(inputData, "NmbUnusedShowers",            axisTitles = "Number of Unused Showers",                       binning = (11, -0.5, 10.5), toposToPlot = toposToPlot, **kwargs)
        # overlayTopologies(inputData, "EnergyUnusedShowers",         axisTitles = "Unused Shower Energy (GeV)",                     binning = (60, 0, 6),       toposToPlot = toposToPlot, **kwargs)
        overlayTopologies(inputData, "BestMissingMatchDistTOF",     axisTitles = "Distance to best ToF match (cm)",                binning = (25, 0, 250),     toposToPlot = toposToPlot, **kwargs)
        # overlayTopologies(inputData, "BestMissingMatchDistBCAL",    axisTitles = "Distance to best BCAL match (cm)",               binning = (20, 0, 200),     toposToPlot = toposToPlot, **kwargs)
        # overlayTopologies(inputData, "MissingMassSquared",          axisTitles = "(#it{m}^{miss}_{kin. fit})^{2} (GeV/#it{c}^{2})^{2}", binning = (125, -0.5, 4.5), toposToPlot = toposToPlot, **kwargs)
        overlayTopologies(inputData, "MissingMassSquared_Measured", axisTitles = "(#it{m}^{miss}_{measured})^{2} (GeV/#it{c}^{2})^{2}", binning = (125, -0.5, 4.5), toposToPlot = toposToPlot, **kwargs)

  # the histograms below are filled for all data types
  cutsArgs = [
    {},  # no extra cut
    {"additionalFilter" : "(NmbUnusedShowers == 0)", "pdfFileNameSuffix" : "_noUnusedShowers"},
  ]
  for kwargs in cutsArgs:

    plot1D(inputData, "AccidWeightFactor",        axisTitles = "RF Weight",                             binning = (1000, -2, 2),    **kwargs, weightVariable = None)
    plot1D(inputData, "KinFitPVal",               axisTitles = "#it{#chi}^{2}_{kim. fit} #it{P}-value", binning = (150, 0, 1),      **kwargs)
    plot1D(inputData, "NmbUnusedShowers",         axisTitles = "Number of Unused Showers",              binning = (11, -0.5, 10.5), **kwargs)
    plot1D(inputData, "BeamEnergy",               axisTitles = "#it{E}_{beam} (GeV)",                   binning = (180, 3, 12),     **kwargs)
    plot1D(inputData, "BestMissingMatchDistTOF",  axisTitles = "Distance to best ToF match (cm)",       binning = (25, 0, 250),     **kwargs)
    plot1D(inputData, "BestMissingMatchDistBCAL", axisTitles = "Distance to best BCAL match (cm)",      binning = (20, 0, 200),     **kwargs)

    sideBandYTitle = "Number of Combos (RF-Sideband)"
    # sideBandArgs: Dict[str, Any] = {
    #   "weightVariable"    : ("AccidWeightFactorSb", "1 - AccidWeightFactor"),
    #   "additionalFilter"  : kwargs.get("additionalFilter", None),
    #   "pdfFileNameSuffix" : "_Sb" + kwargs.get("pdfFileNameSuffix", ""),
    # }
    # plot1D(inputData, ("MissingMass", "sqrt(MissingMassSquared)"), axisTitles = "#it{m}^{miss}_{kin. fit} (GeV/#it{c}^{2})",                   binning = (100, 0, 2), **kwargs)
    # plot1D(inputData, ("MissingMass", "sqrt(MissingMassSquared)"), axisTitles = "#it{m}^{miss}_{kin. fit} (GeV/#it{c}^{2});" + sideBandYTitle, binning = (100, 0, 2), **sideBandArgs)
    # plot1D(inputData, "MissingMassSquared",  axisTitles = "(#it{m}^{miss}_{kin. fit})^{2} (GeV/#it{c}^{2})^{2}",                   binning = (225, -0.5, 4), **kwargs)
    # plot1D(inputData, "MissingMassSquared",  axisTitles = "(#it{m}^{miss}_{kin. fit})^{2} (GeV/#it{c}^{2})^{2};" + sideBandYTitle, binning = (225, -0.5, 4), **sideBandArgs)
    # plot1D(inputData, ("MissingMass_Measured", "sqrt(MissingMassSquared_Measured)"), axisTitles = "#it{m}^{miss}_{measured} (GeV/#it{c}^{2})",                   binning = (100, 0, 2), **kwargs)
    # plot1D(inputData, ("MissingMass_Measured", "sqrt(MissingMassSquared_Measured)"), axisTitles = "#it{m}^{miss}_{measured} (GeV/#it{c}^{2});" + sideBandYTitle, binning = (100, 0, 2), **sideBandArgs)

    # missing-mass squared distributions
    mm2HistDef:         Dict[str, Any] = {"variable" : "MissingMassSquared_Measured", "axisTitles" : "(#it{m}^{miss}_{measured})^{2} (GeV/#it{c}^{2})^{2}",                   "binning" : (125, -0.5, 4.5)}
    mm2HistDefSideBand: Dict[str, Any] = {"variable" : "MissingMassSquared_Measured", "axisTitles" : "(#it{m}^{miss}_{measured})^{2} (GeV/#it{c}^{2})^{2};" + sideBandYTitle, "binning" : (125, -0.5, 4.5), "weightVariable" : ("AccidWeightFactorSb", "1 - AccidWeightFactor")}
    overlayCases(inputData, **mm2HistDef, **kwargs)
    overlayCases(inputData, **mm2HistDefSideBand, pdfFileNameSuffix = f"_Sb" + kwargs.get("pdfFileNameSuffix", ""), additionalFilter = kwargs.get("additionalFilter", None))
    # plot overall distributions for each case
    for case, caseFilter in FILTER_CASES.items():
      caseData = inputData.Filter(caseFilter)
      plot1D(caseData, **mm2HistDef,         pdfFileNameSuffix = f"_{case}"    + kwargs.get("pdfFileNameSuffix", ""), additionalFilter = kwargs.get("additionalFilter", None))
      plot1D(caseData, **mm2HistDefSideBand, pdfFileNameSuffix = f"_{case}_Sb" + kwargs.get("pdfFileNameSuffix", ""), additionalFilter = kwargs.get("additionalFilter", None))
    # kinematicBinnings  = [
    #   # beam energy
    #   # {"variable" : "BeamEnergy",         "label" : "Beam Energy",                   "unit" : "GeV",   "nmbBins" :  9, "range" : (3.0, 12.0)},  # spring 2017
    #   {"variable" : "BeamEnergy",         "label" : "Beam Energy",                   "unit" : "GeV",   "nmbBins" : 10, "range" : (5.5, 11.5)},  # spring 2018
    #   # momentum of missing proton
    #   {"variable" : "MissingProtonP",     "label" : "#it{p}^{miss}_{kin. fit}",      "unit" : "GeV/#it{c}", "nmbBins" : 10, "range" : (0, 3.5)},
    #   # polar angle of missing proton
    #   {"variable" : "MissingProtonTheta", "label" : "#it{#theta}^{miss}_{kin. fit}", "unit" : "deg",   "nmbBins" : 10, "range" : (0, 65)},
    #   # azimuthal angle of missing proton
    #   {"variable" : "MissingProtonPhi",   "label" : "#it{#phi}^{miss}_{kin. fit}",   "unit" : "deg",   "nmbBins" : 10, "range" : (-180, +180)},
    # ]
    # for kinematicBinning in kinematicBinnings:
    #   kinBinVariable = kinematicBinning["variable"]
    #   nmbKinBins     = kinematicBinning["nmbBins"]
    #   kinBinRange    = kinematicBinning["range"]
    #   kinBinWidth    = (kinBinRange[1] - kinBinRange[0]) / float(nmbKinBins)
    #   # plot distributions for kinematic bins
    #   for kinBinIndex in range(nmbKinBins):
    #     kinBinMin = kinBinRange[0] + kinBinIndex * kinBinWidth
    #     kinBinMax = kinBinMin + kinBinWidth
    #     kinBinFilter = f"(({kinBinMin} < {kinBinVariable}) and ({kinBinVariable} < {kinBinMax}))"
    #     kinBinData = inputData.Filter(kinBinFilter)
    #     overlayCases(kinBinData, **mm2HistDef, pdfFileNameSuffix = f"_{kinBinVariable}_{kinBinMin}_{kinBinMax}" + kwargs.get("pdfFileNameSuffix", ""), additionalFilter = kwargs.get("additionalFilter", None))

    plot1D(inputData, "MissingProtonP",     axisTitles = "#it{p}^{miss}_{kin. fit} (GeV/#it{c})", binning = (500, 0, 10),     **kwargs)
    plot1D(inputData, "MissingProtonTheta", axisTitles = "#it{#theta}^{miss}_{kin. fit} (deg)",   binning = (200, 0, 100),    **kwargs)
    plot1D(inputData, "MissingProtonPhi",   axisTitles = "#it{#phi}^{miss}_{kin. fit} (deg)",     binning = (180, -180, 180), **kwargs)

    plot2D(inputData, xVariable = "MissingProtonTheta", yVariable = "MissingProtonP",   axisTitles = "#it{#theta}^{miss}_{kin. fit} (deg);#it{p}^{miss}_{kin. fit} (GeV/#it{c})", binning = (180, 0, 90, 400, 0, 9),      **kwargs)
    plot2D(inputData, xVariable = "MissingProtonTheta", yVariable = "MissingProtonPhi", axisTitles = "#it{#theta}^{miss}_{kin. fit} (deg);#it{#phi}^{miss}_{kin. fit} (deg)",     binning = (180, 0, 90, 180, -180, 180), **kwargs)
    # plot2D(inputData, xVariable = "MissingProtonTheta_Measured", yVariable = "MissingProtonP_Measured",   axisTitles = "#it{#theta}^{miss}_{measured} (deg);#it{p}^{miss}_{measured} (GeV/#it{c})", binning = (180, 0, 90, 400, 0, 9),      **kwargs)
    # plot2D(inputData, xVariable = "MissingProtonTheta_Measured", yVariable = "MissingProtonPhi_Measured", axisTitles = "#it{#theta}^{miss}_{measured} (deg);#it{#phi}^{miss}_{measured} (deg)",     binning = (180, 0, 90, 360, -180, 180), **kwargs)

    # plot1D(inputData, "UnusedDeltaP",      axisTitles = "#it{p}^{miss}_{unused} #minus #it{p}^{miss}_{kin. fit} (GeV/#it{c})",                 binning = (600, -6, 6),     **kwargs)
    plot1D(inputData, "UnusedDeltaPOverP", axisTitles = "(#it{p}^{miss}_{unused} #minus #it{p}^{miss}_{kin. fit}) / #it{p}^{miss}_{kin. fit}", binning = (500, -2, 2),     **kwargs)
    plot1D(inputData, "UnusedDeltaTheta",  axisTitles = "#it{#theta}^{miss}_{unused} #minus #it{#theta}^{miss}_{kin. fit} (deg)",              binning = (200, -100, 100), **kwargs)
    plot1D(inputData, "UnusedDeltaPhi",    axisTitles = "#it{#phi}^{miss}_{unused} #minus #it{#phi}^{miss}_{kin. fit} (deg)",                  binning = (360, -180, 180), **kwargs)
    # overlayCases(inputData, "UnusedDeltaP",      axisTitles = "#it{p}^{miss}_{unused} #minus #it{p}^{miss}_{kin. fit} (GeV/#it{c})",                 binning = (600, -6, 6),     **kwargs)
    overlayCases(inputData, "UnusedDeltaPOverP", axisTitles = "(#it{p}^{miss}_{unused} #minus #it{p}^{miss}_{kin. fit}) / #it{p}^{miss}_{kin. fit}", binning = (500, -2, 2),     **kwargs)
    overlayCases(inputData, "UnusedDeltaTheta",  axisTitles = "#it{#theta}^{miss}_{unused} #minus #it{#theta}^{miss}_{kin. fit} (deg)",              binning = (200, -100, 100), **kwargs)
    overlayCases(inputData, "UnusedDeltaPhi",    axisTitles = "#it{#phi}^{miss}_{unused} #minus #it{#phi}^{miss}_{kin. fit} (deg)",                  binning = (360, -180, 180), **kwargs)
    # unusedTrackData = inputData.Filter("(NmbUnusedTracks == 1)")  # make sure unused track info exists; NOTE! this assumes that there is maximum 1 unused track
    # plot2D(unusedTrackData, xVariable = ("UnusedP_",     "UnusedP[0]"),     yVariable = "MissingProtonP",     axisTitles = "#it{p}^{miss}_{unused} (GeV/#it{c});#it{p}^{miss}_{kin. fit} (GeV/#it{c})", binning = (400, 0, 9, 400, 0, 9),           **kwargs)
    # plot2D(unusedTrackData, xVariable = ("UnusedTheta_", "UnusedTheta[0]"), yVariable = "MissingProtonTheta", axisTitles = "#it{#theta}^{miss}_{unused} (deg);#it{#theta}^{miss}_{kin. fit} (deg)",     binning = (360, 0, 180, 360, 0, 180),       **kwargs)
    # plot2D(unusedTrackData, xVariable = ("UnusedPhi_",   "UnusedPhi[0]"),   yVariable = "MissingProtonPhi",   axisTitles = "#it{#phi}^{miss}_{unused} (deg);#it{#phi}^{miss}_{kin. fit} (deg)",         binning = (360, -180, 180, 360, -180, 180), **kwargs)
