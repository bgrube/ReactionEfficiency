#!/usr/bin/env python3


from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
import functools
import numpy as np
from typing import Any

import ROOT

from plotFitResults import (
  getAxisInfoForBinningVar,
  plotGraphs1D,
)
from plotTools import (
  callMemberFunctionsWithArgs,
  drawZeroLine,
  getCbFriendlyRootColor,
  makeDirPath,
  printGitInfo,
  setupPlotStyle,
)


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
FILTER_CASES: dict[str, str] = {
  "Total"   : "(true)",
  "Found"   : "(TrackFound == true)",
  "Missing" : "(TrackFound == false)",
}

COLOR_CASES = {
  "Total"   : ROOT.kGray,
  "Found"   : ROOT.kGreen + 2,
  "Missing" : ROOT.kRed + 1,
}


def overlayDataSamples1D(
  dataSamples:       dict[str, dict[str, Any]],  # TFile or RDataFrame and style definitions for each data-set label
  histName:          str | None                      = None,  # name of histogram to plot; required for TFile
  variable:          str | tuple[str, str] | None    = None,        # variable to plot; may be column name, or tuple with new column definition; required for RDataFrame
  axisTitles:        str | None                      = None,  # semicolon-separated list; required for RDataFrame
  binning:           tuple[int, float, float] | None = None,  # tuple with binning definition; required for RDataFrame
  weightVariable:    str | tuple[str, str] | None    = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
  pdfFileNamePrefix: str                             = "Proton_4pi_",
  pdfFileNameSuffix: str                             = "",
  pdfDirName:        str                             = "./",
  additionalFilter:  str | None                      = None,
  histTitle:         str | None                      = None,
) -> None:
  """Overlays 1D histograms generated from given trees or read from given files"""
  print("Overlaying " + (f"histograms '{histName}'" if histName else f"distributions for '{variable}'") + f" for data samples '{', '.join(dataSamples.keys())}'")
  pdfFileBaseName = histName if histName is not None else variable
  hStack = ROOT.THStack(pdfFileBaseName, ("" if histTitle is None else histTitle) + f";{setDefaultYAxisTitle(axisTitles)}")
  hists: list[ROOT.TH1D] = []  # keep histograms in memory
  normIntegral = None  # index of histogram to normalize to
  for dataLabel, dataSample in dataSamples.items():
    hist: ROOT.TH1D | None = None
    if "TFile" in dataSample:
      # read histogram from file
      assert histName is not None, f"Name of histogram to read from file '{dataSample['TFile'].GetPath()}' required."
      hist = dataSample["TFile"].Get(histName)
      hist.SetTitle(dataLabel)
    elif "RDataFrame" in dataSample:
      assert variable is not None and axisTitles is not None and binning is not None, f"Need variable name (={variable}), axis titles (={axisTitles}), and binning (={binning})."
      hist = getHistND(dataSample["RDataFrame"], (variable,), setDefaultYAxisTitle(axisTitles), binning, weightVariable,
                       filterExpression = additionalFilter, histNameSuffix = dataLabel, histTitle = dataLabel).GetPtr()
    else:
      raise KeyError(f"Data sample must contain either 'TFile' or 'RDataFrame' key: {dataSample}")
    assert hist is not None, "Could not create histogram"
    callMemberFunctionsWithArgs(hist, dataSample)
    if dataSample.get("normToThis", False):
      print(f"Normalizing all histograms to '{dataLabel}'")
      normIntegral = hist.Integral()
    hists.append(hist)
    hStack.Add(hist)
  # normalize histograms
  if normIntegral is not None:
    for hist in hists:
      hist.Scale(normIntegral / hist.Integral())
  # draw distributions
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{pdfFileBaseName}_overlay_{'_'.join(dataSamples.keys())}{pdfFileNameSuffix}")
  hStack.Draw("NOSTACK HIST")
  # add legend
  # canv.BuildLegend()  # automatic placement with width 0.3 and height 0.21
  canv.BuildLegend(0.3, 0.15, 0.3, 0.15)  # automatic placement with width 0.3 and height 0.15
  # canv.BuildLegend(0.7, 0.65, 0.99, 0.99)
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


def drawHistogram(
  inFileName:        str,
  histName:          str,
  rebinFactor:       int | Sequence[int] = 1,  # if integer -> rebin x axis; if sequence of 2 integers -> rebin x and y axes
  drawOption:        str                 = "HIST",
  pdfFileNamePrefix: str                 = "Proton_4pi_",
  pdfFileNameSuffix: str                 = "",
  pdfDirName:        str                 = "./",
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
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


def getHistND(
  inputData:        ROOT.RDataFrame,
  variables:        tuple[str | tuple[str, str], ...],  # variable(s) to plot; is tuple of either column names or tuples with new column definitions; defines dimension of histogram
  axisTitles:       str,                                # semicolon-separated list
  binning:          tuple[int, float, float] | tuple[int, float, float, int, float, float],  # tuple with 1D or 2D binning definitions
  weightVariable:   str | tuple[str, str] | None = None,  # may be None (= no weighting), string with column name, or tuple with new column definition
  filterExpression: str | None                   = None,
  histNameSuffix:   str                          = "",
  histTitle:        str                          = "",
) -> ROOT.TH1D | ROOT.TH2D:
  """Creates histogram from given variables in RDataFrame, applying optional weighting and filtering"""
  histDim = len(variables)
  assert 1 <= histDim <= 2, "currently, only 1D and 2D histograms are supported"
  # apply additional filters, if defined
  data = inputData.Filter(filterExpression) if filterExpression else inputData
  columnNames: list[str] = [""] * len(variables)
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
  axisTitles:    str | None,  # semicolon-separated list
  defaultYTitle: str = "Number of Combos (RF-subtracted)",
) -> str:
  """Sets default y-axis title if not provided by `axisTitles`"""
  if (axisTitles is None):
    return ";" + defaultYTitle
  titles = axisTitles.split(";")
  if (len(titles) == 1):
    return titles[0] + ";" + defaultYTitle
  elif (len(titles) == 0):
    return ";" + defaultYTitle
  else:
    return axisTitles


def plot1D(
  inputData:         ROOT.RDataFrame,
  variable:          str | tuple[str, str],     # variable to plot; may be column name, or tuple with new column definition
  axisTitles:        str,                       # semicolon-separated list
  binning:           tuple[int, float, float],  # tuple with binning definition
  weightVariable:    str | tuple[str, str] | None = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
  pdfFileNamePrefix: str                          = "Proton_4pi_",
  pdfFileNameSuffix: str                          = "",
  pdfDirName:        str                          = "./",
  additionalFilter:  str | None                   = None,
) -> None:
  """Plots 1D distribution for given variable, applying optional weighting and filtering"""
  hist: ROOT.TH1D = getHistND(inputData, (variable,), setDefaultYAxisTitle(axisTitles), binning, weightVariable, filterExpression = additionalFilter)
  # draw distributions
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{hist.GetName()}{pdfFileNameSuffix}")
  hist.Draw("HIST")
  drawZeroLine(hist)
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


def plot2D(
  inputData:         ROOT.RDataFrame,
  xVariable:         str | tuple[str, str],  # x variable to plot; may be column name, or tuple with new column definition
  yVariable:         str | tuple[str, str],  # y variable to plot; may be column name, or tuple with new column definition
  axisTitles:        str,                    # semicolon-separated list
  binning:           tuple[int, float, float, int, float, float],  # tuple with binning definition
  weightVariable:    str | tuple[str, str] | None = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
  pdfFileNamePrefix: str                          = "Proton_4pi_",
  pdfFileNameSuffix: str                          = "",
  pdfDirName:        str                          = "./",
  additionalFilter:  str | None                   = None,
) -> None:
  """Plots 2D distribution for given x and y variables, applying optional weighting and filtering"""
  hist: ROOT.TH2D = getHistND(inputData, (xVariable, yVariable), axisTitles, binning, weightVariable, filterExpression = additionalFilter)
  # draw distributions
  canv = ROOT.TCanvas(f"{pdfFileNamePrefix}{hist.GetName()}{pdfFileNameSuffix}")
  hist.Draw("COLZ")
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


def getResolutionGraph(
  inputData:           ROOT.RDataFrame,
  diffVariable:        str,                       # residual from which resolution is measured
  diffVariableBinning: tuple[int, float, float],  # tuple with binning definition for diffVariable
  resBinningVariable:  str | tuple[str, str],     # variable to bin resolution in; may be column name, or tuple with new column definition
  resBinning:          tuple[int, float, float],  # tuple with binning definition for resBinningVariable
  weightVariable:      str | tuple[str, str] | None = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
  additionalFilter:    str | None                   = None,
  histTitle:           str | None                   = None,
  histPdfFileName:     str | None                   = None,
) -> ROOT.TGraphErrors:
  """Plots resolution of given variable in the in bins of `resBinningVariable`"""
  # !Note! resolution variables are arrays with length NmbTruthTracks, need to define dummy column for first element
  hist: ROOT.TH2D = getHistND(inputData, ((f"{diffVariable}_", f"{diffVariable}[0]"), resBinningVariable), axisTitles = "",
                              binning = (*diffVariableBinning, *resBinning), weightVariable = weightVariable, filterExpression = additionalFilter)
  # draw distributions
  if histPdfFileName is not None:
    canv = ROOT.TCanvas(histPdfFileName)
    canv.SetLogz()
    if histTitle is not None:
      hist.SetTitle(histTitle)
    hist.Draw("COLZ")
    canv.SaveAs(histPdfFileName)
  # construct resolution graph
  resBinningAxis = hist.GetYaxis()
  xVals = np.array([resBinningAxis.GetBinCenter(resBinIndex)       for resBinIndex in range(1, resBinningAxis.GetNbins() + 1)], dtype = np.float64)
  xErrs = np.array([resBinningAxis.GetBinWidth (resBinIndex) / 2.0 for resBinIndex in range(1, resBinningAxis.GetNbins() + 1)], dtype = np.float64)
  yVals = np.array([], dtype = np.float64)
  yErrs = np.array([], dtype = np.float64)
  for resBinIndex in range(1, resBinningAxis.GetNbins() + 1):
    proj = hist.ProjectionX("_px", resBinIndex, resBinIndex, "E")
    yVals = np.append(yVals, (proj.GetStdDev(),))
    yErrs = np.append(yErrs, (proj.GetStdDevError(),))
  return ROOT.TGraphErrors(len(xVals), xVals, yVals, xErrs, yErrs)


def overlayResolutions(
  inputData:             Sequence[tuple[str, ROOT.RDataFrame]],
  resVariableName:       str,                       # name of variable for which resolution is estimated
  diffVariable:          str,                       # residual from which resolution is measured
  diffVariableBinning:   tuple[int, float, float],  # tuple with binning definition for diffVariable
  resBinningVariable:    str | tuple[str, str],     # variable to bin resolution in; may be column name, or tuple with new column definition
  resBinning:            tuple[int, float, float],  # tuple with binning definition for resBinningVariable
  weightVariable:        str | tuple[str, str] | None      = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
  pdfFileNamePrefix:     str                               = "Proton_4pi_",
  pdfFileNameSuffix:     str                               = "",
  pdfDirName:            str                               = "./",
  diffVariableAxisTitle: str | None                        = None,
  additionalFilter:      str | None                        = None,
  resPlotRange:          tuple[float | None, float | None] = (None, None),
) -> None:
  """Plots resolution of given variable in the in bins of `resBinningVariable`"""
  resBinningVariableName, resBinningVariableLabel, resBinningVariableUnit = getAxisInfoForBinningVar(resBinningVariable)
  resGraphs: list[tuple[str, ROOT.TGraphErrors]] = []
  for label, data in inputData:
    resGraphs.append((label,
      getResolutionGraph(
        inputData           = data,
        diffVariable        = diffVariable,
        diffVariableBinning = diffVariableBinning,
        resBinningVariable  = resBinningVariable,
        resBinning          = resBinning,
        weightVariable      = weightVariable,
        additionalFilter    = additionalFilter,
        histTitle           = f";{'' if diffVariableAxisTitle is None else diffVariableAxisTitle};{resBinningVariableLabel} ({resBinningVariableUnit})",
        histPdfFileName     = f"{pdfDirName}/{pdfFileNamePrefix}resolution2D_{diffVariable}_{resBinningVariable}_{label}{pdfFileNameSuffix}.pdf",
      )))
  _, resVariableLabel, resVariableUnit = getAxisInfoForBinningVar(resVariableName)
  plotGraphs1D(
    graphOrGraphs     = resGraphs,
    binningVar        = resBinningVariableName,
    yAxisTitle        = f"#it{{#sigma}}_{{{resVariableLabel}}} ({resVariableUnit})",
    pdfDirName        = pdfDirName,
    pdfFileBaseName   = f"resolution_{resVariableName}",
    pdfFileNamePrefix = pdfFileNamePrefix,
    pdfFileNameSuffix = pdfFileNameSuffix,
    graphMinimum      = resPlotRange[0],
    graphMaximum      = resPlotRange[1],
    skipBlack         = True if len(resGraphs) > 1 else False,
  )


def overlayCases(
  inputData:         ROOT.RDataFrame,
  variable:          str | tuple[str, str],     # variable to plot; may be column name, or tuple with new column definition
  axisTitles:        str,                       # semicolon-separated list
  binning:           tuple[int, float, float],  # tuple with binning definition
  weightVariable:    str | tuple[str, str] | None = "AccidWeightFactor",  # may be None (= no weighting), string with column name, or tuple with new column definition
  pdfFileNamePrefix: str                          = "Proton_4pi_",
  pdfFileNameSuffix: str                          = "",
  pdfDirName:        str                          = "./",
  additionalFilter:  str | None                   = None,
) -> None:
  """Overlays 1D distributions of given variable for "Total", "Found", and "Missing" cases"""
  data = inputData.Filter(additionalFilter) if additionalFilter else inputData
  # overlay distributions for cases
  hStack = ROOT.THStack(f"{variable}", ";" + setDefaultYAxisTitle(axisTitles))
  hists = []
  for case, caseFilter in FILTER_CASES.items():
    hist: ROOT.TH1D = getHistND(data, (variable,), setDefaultYAxisTitle(axisTitles), binning, weightVariable,
                                filterExpression = caseFilter, histNameSuffix = case, histTitle = case)
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
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


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
  weightVariable:   str | tuple[str, str] | None = None,  # may be None (= no weighting), string with column name, or tuple with new column definition
  filterExpression: str | None                   = None,
  histNameSuffix:   str = "",
) -> tuple[list[str], ROOT.TH1F]:
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


def getCategoricalTH1AsDict(hist: ROOT.TH1) -> dict[str, float]:
  """Returns categorical histogram as dict { bin label : bin content }"""
  xAxis = hist.GetXaxis()
  return {xAxis.GetBinLabel(binIndex) : hist.GetBinContent(binIndex)
    for binIndex in range(1, xAxis.GetNbins() + 1)}


def plotTopologyHist(
  inputData:         ROOT.RDataFrame,
  normalize:         bool       = False,
  maxNmbTopologies:  int        = 10,
  additionalFilter:  str | None = None,
  pdfFileNamePrefix: str        = "Proton_4pi_",
  pdfFileNameSuffix: str        = "",
  pdfDirName:        str        = "./",
) -> None:
  """Plots categorical histogram with counts or fraction for each generated topology, applying optional weighting and filtering"""
  # get histogram data
  topoNames: dict[str, list[str]] = {}  # dictionary of ordered list of topology names { case : [ topologyName ] }
  topoHists: dict[str, ROOT.TH1F] = {}  # dictionary of topology histograms { case : topologyHist }
  for case, caseFilter in FILTER_CASES.items():
    caseData = inputData.Filter(caseFilter)
    topoNames[case], topoHists[case] = getTopologyHist(caseData, weightVariable = "AccidWeightFactor", filterExpression = additionalFilter, histNameSuffix = case + ("_norm" if normalize else ""))
  # overlay distributions for cases
  hStack = ROOT.THStack(f"topologies",  ";;" + ("Fraction" if normalize else "Number") + " of Combos (RF-subtracted)" + (" [%]" if normalize else ""))
  topoLabels = topoNames["Total"]
  hists: dict[str, ROOT.TH1] = {}  # memorize plots to print
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
    integral = hist.Integral()
    if normalize and integral != 0:
      hist.Scale(100 / integral)
    hist.SetLineColor(COLOR_CASES[case])
    if case == "Total":
      hist.SetFillColor(COLOR_CASES[case])
    hists[case] = hist
    hStack.Add(hist)
    print(f"plotTopologyHist(): {case} " + (f"purity = {hist.GetBinContent(1)}%" if normalize else f"signal = {hist.GetBinContent(1)} of {integral} combos"))
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
  canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


def overlayTopologies(
  inputData:         ROOT.RDataFrame,
  variable:          str | tuple[str, str],     # variable to plot; may be column name, or tuple with new column definition
  axisTitles:        str,                       # semicolon-separated list
  binning:           tuple[int, float, float],  # tuple with binning definition
  toposToPlot:       dict[str, list[str]],      # topologies to plot for each case
  additionalFilter:  str | None = None,
  pdfFileNamePrefix: str        = "Proton_4pi_",
  pdfFileNameSuffix: str        = "_MCbggen_topologies",
  pdfDirName:        str        = "./",
) -> None:
  """Overlays 1D distributions for given variable from overall data sample and distributions for the `maxNmbTopologies` topologies with the largest number of combos from the bggen MC sample"""
  data = inputData.Filter(additionalFilter) if additionalFilter else inputData
  for case, caseFilter in FILTER_CASES.items():
    caseData = data.Filter(caseFilter)
    # get topologies with largest number of combos for given case
    hStack = ROOT.THStack(f"{variable}_{case}", f"{case};{setDefaultYAxisTitle(axisTitles)}")
    hists = []
    # overlay distributions for topologies
    for index, topo in enumerate(toposToPlot[case]):
      hist: ROOT.TH1D = getHistND(caseData, (variable,), setDefaultYAxisTitle(axisTitles), binning, "AccidWeightFactor",
                                  filterExpression = (f'ThrownTopology.GetString() == "{topo}"' if topo != "Total" else "true"), histNameSuffix = f"{case}_{topo}", histTitle = topo)
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
    canv.SaveAs(f"{pdfDirName}/{canv.GetName()}.pdf")


def makeKinematicPlotsOverlays(
  dataSamples: dict[str, dict[str, Any]],  # RDataFrame and style definitions for each data-set label
  pdfDirName:  str = "./",
) -> None:
  """Overlays kinematic distributions of given data samples"""
  # overlayDataSamples1D(dataSamples, variable = "NmbUnusedShowers", axisTitles = "Number of Unused Showers", binning = (11, -0.5, 10.5), pdfDirName = pdfDirName)
  kwargss = (
    {
      "additionalFilter"  : "(NmbUnusedShowers == 0)",
      "pdfFileNameSuffix" : "_noUnusedShowers",
      "pdfDirName"        : pdfDirName,
    },
    {
      "additionalFilter"  : "(NmbUnusedShowers == 0) && ((5.5 < BeamEnergy) && (BeamEnergy < 11.0))",
      "pdfFileNameSuffix" : "_noUnusedShowers_beamEnergy_5.5_11.0",
      "pdfDirName"        : pdfDirName,
    },
    {
      "additionalFilter"  : "(NmbUnusedShowers == 0) && ((7.5 < BeamEnergy) && (BeamEnergy < 9.0))",
      "pdfFileNameSuffix" : "_noUnusedShowers_beamEnergy_7.5_9.0",
      "pdfDirName"        : pdfDirName,
    },
  )
  for kwargs in kwargss:
    overlayDataSamples1D(dataSamples, variable = "BeamEnergy",         axisTitles = "#it{E}_{beam} (GeV)",                   binning = (180,    3,  12), **kwargs)
    overlayDataSamples1D(dataSamples, variable = "KinFitPVal",         axisTitles = "#it{#chi}^{2}_{kin. fit} #it{P}-value", binning = (150,    0,   1), **kwargs)
    overlayDataSamples1D(dataSamples, variable = "MissingProtonP",     axisTitles = "#it{p}_{miss}^{kin. fit} (GeV/#it{c})", binning = (250,    0,   5), **kwargs)
    overlayDataSamples1D(dataSamples, variable = "MissingProtonTheta", axisTitles = "#it{#theta}_{miss}^{kin. fit} (deg)",   binning = (200,    0, 100), **kwargs)
    overlayDataSamples1D(dataSamples, variable = "MissingProtonPhi",   axisTitles = "#it{#phi}_{miss}^{kin. fit} (deg)",     binning = (180, -180, 180), **kwargs)
    overlayDataSamples1D(dataSamples, variable = "FourPiMass",         axisTitles = "#it{m}_{#it{#pi}^{#plus}#it{#pi}^{#minus}#it{#pi}^{#plus}#it{#pi}^{#minus}} (GeV/#it{c}^{2})", binning = (200, 0, 5), **kwargs)
  # unused track
  overlayDataSamples1D(dataSamples, variable = "UnusedDeltaPOverP", axisTitles = "(#it{p}_{miss}^{unused} #minus #it{p}_{miss}^{kin. fit}) / #it{p}_{miss}^{kin. fit}", binning = (375, -1.5, +1.5), **kwargss[0])
  overlayDataSamples1D(dataSamples, variable = "UnusedDeltaTheta",  axisTitles = "#it{#theta}_{miss}^{unused} #minus #it{#theta}_{miss}^{kin. fit} (deg)",              binning = (100, -50,  +50),  **kwargss[0])
  overlayDataSamples1D(dataSamples, variable = "UnusedDeltaPhi",    axisTitles = "#it{#phi}_{miss}^{unused} #minus #it{#phi}_{miss}^{kin. fit} (deg)",                  binning = (200, -100, +100), **kwargss[0])
  # missing mass squared
  for case, caseFilter in FILTER_CASES.items():
    kwargs = {
      "additionalFilter"  : f"((NmbUnusedShowers == 0) and {caseFilter})",
      "pdfFileNameSuffix" : f"_{case}_noUnusedShowers",
      "pdfDirName"        : pdfDirName,
    }
    overlayDataSamples1D(dataSamples, variable = "MissingMassSquared_Measured",
                          axisTitles = "(#it{m}_{miss}^{meas.})^{2} (GeV/#it{c}^{2})^{2}", binning = (125, -0.5, 4.5), **kwargs)


def makeKinematicPlotsMc(
  dataSample: ROOT.RDataFrame,
  isMcBggen:  bool,
  pdfDirName: str = "./",
) -> None:
  """Plots kinematic distributions for given Monte Carlo data"""
  filterTopologies = {
    ""                                             : None,
    "__sig"                                        : '(ThrownTopology.GetString() == "2#pi^{#plus}2#pi^{#minus}p")',
    # "__2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0}]" : '(ThrownTopology.GetString() == "2#gamma2#pi^{#plus}2#pi^{#minus}p[#pi^{0}]")',
    "__bkg"                                        : '(ThrownTopology.GetString() != "2#pi^{#plus}2#pi^{#minus}p")',
  }
  for suffix, filter in filterTopologies.items():
    kwargs = {
      "additionalFilter"  : "(NmbTruthTracks == 1)" + ("" if filter is None else f" && {filter}"),
      "pdfFileNameSuffix" : suffix,
      "pdfDirName"        : pdfDirName,
    }
    overlayCases(dataSample, "TruthDeltaP",      axisTitles = "#it{p}_{miss}^{truth} #minus #it{p}_{miss}^{kin. fit} (GeV/#it{c})",                 binning = (600,   -4,   +4), **kwargs)
    overlayCases(dataSample, "TruthDeltaPOverP", axisTitles = "(#it{p}_{miss}^{truth} #minus #it{p}_{miss}^{kin. fit}) / #it{p}_{miss}^{kin. fit}", binning = (500,   -2,   +2), **kwargs)
    overlayCases(dataSample, "TruthDeltaTheta",  axisTitles = "#it{#theta}_{miss}^{truth} #minus #it{#theta}_{miss}^{kin. fit} (deg)",              binning = (200,  -60,  +60), **kwargs)
    overlayCases(dataSample, "TruthDeltaPhi",    axisTitles = "#it{#phi}_{miss}^{truth} #minus #it{#phi}_{miss}^{kin. fit} (deg)",                  binning = (360, -180, +180), **kwargs)

  if isMcBggen:  # fill plots for bggen Monte Carlo

    cutsArgs: list[dict[str, Any]] = [
      {},  # no extra cut
      {"additionalFilter" : "(NmbUnusedShowers == 0)",                                      "pdfFileNameSuffix" : f"_noUnusedShowers"},  # no unused showers
      # {"additionalFilter" : "((NmbUnusedShowers == 0) and (BestMissingMatchDistTOF < 40))", "pdfFileNameSuffix" : f"_noUnusedShowersMatchToF"},  # no unused showers and ToF hit within certain distance
    ]
    for kwargs in cutsArgs:
      plotTopologyHist(dataSample, normalize = False, pdfDirName = pdfDirName, **kwargs)
      plotTopologyHist(dataSample, normalize = True,  pdfDirName = pdfDirName, **kwargs)

    cutsArgs = [
      {},  # no extra cut
      {"additionalFilter" : "(NmbUnusedShowers == 0)",                                      "pdfFileNameSuffix" : "_noUnusedShowers"},  # no unused showers
      # # the two cuts below are equivalent to the one above
      # {"additionalFilter" : "(EnergyUnusedShowers == 0)",                                   "pdfFileNameSuffix" : "_noEnergyUnusedShowers"},
      # {"additionalFilter" : "(NmbUnusedShowers == 0) and (EnergyUnusedShowers == 0)",       "pdfFileNameSuffix" : "_noShowers"}
      # {"additionalFilter" : "((NmbUnusedShowers == 0) and (BestMissingMatchDistTOF < 40))", "pdfFileNameSuffix" : "_noUnusedShowersMatchToF"},  # no unused showers and ToF hit within certain distance
    ]
    for kwargs in cutsArgs:
      kwargs.update({"pdfDirName" : pdfDirName})
      # get topologies with the largest number of combos for given case
      toposToPlot: dict[str, list[str]] = {}
      for case, caseFilter in FILTER_CASES.items():
        caseData = dataSample.Filter(caseFilter)
        toposToPlot[case], _ = getTopologyHist(caseData, filterExpression = kwargs.get("additionalFilter", None))
        toposToPlot[case] = ["Total"] + toposToPlot[case][:maxNmbTopologies]
      overlayTopologies(dataSample, "NmbUnusedShowers",            axisTitles = "Number of Unused Showers",                       binning = (11, -0.5, 10.5), toposToPlot = toposToPlot, **kwargs)
      # overlayTopologies(dataSample, "EnergyUnusedShowers",         axisTitles = "Unused Shower Energy (GeV)",                     binning = (60, 0, 6),       toposToPlot = toposToPlot, **kwargs)
      # overlayTopologies(dataSample, "BestMissingMatchDistTOF",     axisTitles = "Distance to best ToF match (cm)",                binning = (25, 0, 250),     toposToPlot = toposToPlot, **kwargs)
      # overlayTopologies(dataSample, "BestMissingMatchDistBCAL",    axisTitles = "Distance to best BCAL match (cm)",               binning = (20, 0, 200),     toposToPlot = toposToPlot, **kwargs)
      # overlayTopologies(dataSample, "MissingMassSquared",          axisTitles = "(#it{m}_{miss}^{kin. fit})^{2} (GeV/#it{c}^{2})^{2}", binning = (125, -0.5, 4.5), toposToPlot = toposToPlot, **kwargs)
      overlayTopologies(dataSample, "MissingMassSquared_Measured", axisTitles = "(#it{m}_{miss}^{meas.})^{2} (GeV/#it{c}^{2})^{2}", binning = (125, -0.5, 4.5), toposToPlot = toposToPlot, **kwargs)


def makeKinematicPlotsData(
  dataSample: ROOT.RDataFrame,
  pdfDirName: str = "./",
) -> None:
  """Plots kinematic distributions for given Monte Carlo data"""
  cutsArgs: list[dict[str, Any]] = [
    {},  # no extra cut
    {"additionalFilter" : "(NmbUnusedShowers == 0)", "pdfFileNameSuffix" : "_noUnusedShowers"},
  ]
  for kwargs in cutsArgs:
    kwargs.update({"pdfDirName" : pdfDirName})

    plot1D(dataSample, "AccidWeightFactor",        axisTitles = "RF Weight",                             binning = (1000, -2, 2),    **kwargs, weightVariable = None)
    plot1D(dataSample, "KinFitPVal",               axisTitles = "#it{#chi}^{2}_{kin. fit} #it{P}-value", binning = (150, 0, 1),      **kwargs)
    plot1D(dataSample, "NmbUnusedShowers",         axisTitles = "Number of Unused Showers",              binning = (11, -0.5, 10.5), **kwargs)
    plot1D(dataSample, "BeamEnergy",               axisTitles = "#it{E}_{beam} (GeV)",                   binning = (180, 3, 12),     **kwargs)
    plot1D(dataSample, "BestMissingMatchDistTOF",  axisTitles = "Distance to best ToF match (cm)",       binning = (25, 0, 250),     **kwargs)
    plot1D(dataSample, "BestMissingMatchDistBCAL", axisTitles = "Distance to best BCAL match (cm)",      binning = (20, 0, 200),     **kwargs)
    plot1D(dataSample, "FourPiMass",               axisTitles = "#it{m}_{#it{#pi}^{#plus}#it{#pi}^{#minus}#it{#pi}^{#plus}#it{#pi}^{#minus}} (GeV/#it{c}^{2})", binning = (200, 0, 5), **kwargs)

    sideBandYTitle = "Number of Combos (RF-Sideband)"
    # sideBandArgs: dict[str, Any] = {
    #   "weightVariable"    : ("AccidWeightFactorSb", "1 - AccidWeightFactor"),
    #   "additionalFilter"  : kwargs.get("additionalFilter", None),
    #   "pdfFileNameSuffix" : "_Sb" + kwargs.get("pdfFileNameSuffix", ""),
    #   "pdfDirName"        : pdfDirName,
    # }
    # plot1D(dataSample, ("MissingMass", "sqrt(MissingMassSquared)"), axisTitles = "#it{m}_{miss}^{kin. fit} (GeV/#it{c}^{2})",                   binning = (100, 0, 2), **kwargs)
    # plot1D(dataSample, ("MissingMass", "sqrt(MissingMassSquared)"), axisTitles = "#it{m}_{miss}^{kin. fit} (GeV/#it{c}^{2});" + sideBandYTitle, binning = (100, 0, 2), **sideBandArgs)
    # plot1D(dataSample, "MissingMassSquared",  axisTitles = "(#it{m}_{miss}^{kin. fit})^{2} (GeV/#it{c}^{2})^{2}",                   binning = (225, -0.5, 4), **kwargs)
    # plot1D(dataSample, "MissingMassSquared",  axisTitles = "(#it{m}_{miss}^{kin. fit})^{2} (GeV/#it{c}^{2})^{2};" + sideBandYTitle, binning = (225, -0.5, 4), **sideBandArgs)
    # plot1D(dataSample, ("MissingMass_Measured", "sqrt(MissingMassSquared_Measured)"), axisTitles = "#it{m}_{miss}^{meas.} (GeV/#it{c}^{2})",                   binning = (100, 0, 2), **kwargs)
    # plot1D(dataSample, ("MissingMass_Measured", "sqrt(MissingMassSquared_Measured)"), axisTitles = "#it{m}_{miss}^{meas.} (GeV/#it{c}^{2});" + sideBandYTitle, binning = (100, 0, 2), **sideBandArgs)

    # missing-mass squared distributions
    mm2HistDef:         dict[str, Any] = {"variable" : "MissingMassSquared_Measured", "axisTitles" : "(#it{m}_{miss}^{meas.})^{2} (GeV/#it{c}^{2})^{2}",                   "binning" : (125, -0.5, 4.5)}
    mm2HistDefSideBand: dict[str, Any] = {"variable" : "MissingMassSquared_Measured", "axisTitles" : "(#it{m}_{miss}^{meas.})^{2} (GeV/#it{c}^{2})^{2};" + sideBandYTitle, "binning" : (125, -0.5, 4.5), "weightVariable" : ("AccidWeightFactorSb", "1 - AccidWeightFactor")}
    overlayCases(dataSample, **mm2HistDef, **kwargs)
    overlayCases(dataSample, **mm2HistDefSideBand, pdfFileNameSuffix = f"_Sb" + kwargs.get("pdfFileNameSuffix", ""), additionalFilter = kwargs.get("additionalFilter", None), pdfDirName = pdfDirName)
    # plot overall distributions for each case
    for case, caseFilter in FILTER_CASES.items():
      caseData = dataSample.Filter(caseFilter)
      plot1D(caseData, **mm2HistDef,         pdfFileNameSuffix = f"_{case}"    + kwargs.get("pdfFileNameSuffix", ""), additionalFilter = kwargs.get("additionalFilter", None), pdfDirName = pdfDirName)
      plot1D(caseData, **mm2HistDefSideBand, pdfFileNameSuffix = f"_{case}_Sb" + kwargs.get("pdfFileNameSuffix", ""), additionalFilter = kwargs.get("additionalFilter", None), pdfDirName = pdfDirName)
    # kinematicBinnings  = [
    #   # beam energy
    #   # {"variable" : "BeamEnergy",         "label" : "Beam Energy",                   "unit" : "GeV",   "nmbBins" :  9, "range" : (3.0, 12.0)},  # spring 2017
    #   {"variable" : "BeamEnergy",         "label" : "Beam Energy",                   "unit" : "GeV",   "nmbBins" : 10, "range" : (5.5, 11.5)},  # spring 2018
    #   # momentum of missing proton
    #   {"variable" : "MissingProtonP",     "label" : "#it{p}_{miss}^{kin. fit}",      "unit" : "GeV/#it{c}", "nmbBins" : 10, "range" : (0, 3.5)},
    #   # polar angle of missing proton
    #   {"variable" : "MissingProtonTheta", "label" : "#it{#theta}_{miss}^{kin. fit}", "unit" : "deg",   "nmbBins" : 10, "range" : (0, 65)},
    #   # azimuthal angle of missing proton
    #   {"variable" : "MissingProtonPhi",   "label" : "#it{#phi}_{miss}^{kin. fit}",   "unit" : "deg",   "nmbBins" : 10, "range" : (-180, +180)},
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
    #     kinBinData = dataSample.Filter(kinBinFilter)
    #     overlayCases(kinBinData, **mm2HistDef, pdfFileNameSuffix = f"_{kinBinVariable}_{kinBinMin}_{kinBinMax}" + kwargs.get("pdfFileNameSuffix", ""), additionalFilter = kwargs.get("additionalFilter", None), pdfDirName = pdfDirName)

    plot1D(dataSample, "MissingProtonP",     axisTitles = "#it{p}_{miss}^{kin. fit} (GeV/#it{c})", binning = (250, 0, 5),      **kwargs)
    plot1D(dataSample, "MissingProtonTheta", axisTitles = "#it{#theta}_{miss}^{kin. fit} (deg)",   binning = (200, 0, 100),    **kwargs)
    plot1D(dataSample, "MissingProtonPhi",   axisTitles = "#it{#phi}_{miss}^{kin. fit} (deg)",     binning = (180, -180, 180), **kwargs)

    plot2D(dataSample, xVariable = "MissingProtonTheta", yVariable = "MissingProtonP",   axisTitles = "#it{#theta}_{miss}^{kin. fit} (deg);#it{p}_{miss}^{kin. fit} (GeV/#it{c})", binning = (180, 0, 90, 400, 0, 9),      **kwargs)
    plot2D(dataSample, xVariable = "MissingProtonTheta", yVariable = "MissingProtonPhi", axisTitles = "#it{#theta}_{miss}^{kin. fit} (deg);#it{#phi}_{miss}^{kin. fit} (deg)",     binning = (180, 0, 90, 180, -180, 180), **kwargs)
    # plot2D(dataSample, xVariable = "MissingProtonTheta_Measured", yVariable = "MissingProtonP_Measured",   axisTitles = "#it{#theta}_{miss}^{meas.} (deg);#it{p}_{miss}^{meas.} (GeV/#it{c})", binning = (180, 0, 90, 400, 0, 9),      **kwargs)
    # plot2D(dataSample, xVariable = "MissingProtonTheta_Measured", yVariable = "MissingProtonPhi_Measured", axisTitles = "#it{#theta}_{miss}^{meas.} (deg);#it{#phi}_{miss}^{meas.} (deg)",     binning = (180, 0, 90, 360, -180, 180), **kwargs)

    # plot1D(dataSample, "UnusedDeltaP",      axisTitles = "#it{p}_{miss}^{unused} #minus #it{p}_{miss}^{kin. fit} (GeV/#it{c})",                 binning = (600, -6, 6),     **kwargs)
    plot1D(dataSample, "UnusedDeltaPOverP", axisTitles = "(#it{p}_{miss}^{unused} #minus #it{p}_{miss}^{kin. fit}) / #it{p}_{miss}^{kin. fit}", binning = (500, -2, 2),     **kwargs)
    plot1D(dataSample, "UnusedDeltaTheta",  axisTitles = "#it{#theta}_{miss}^{unused} #minus #it{#theta}_{miss}^{kin. fit} (deg)",              binning = (200, -100, 100), **kwargs)
    plot1D(dataSample, "UnusedDeltaPhi",    axisTitles = "#it{#phi}_{miss}^{unused} #minus #it{#phi}_{miss}^{kin. fit} (deg)",                  binning = (360, -180, 180), **kwargs)
    # overlayCases(dataSample, "UnusedDeltaP",      axisTitles = "#it{p}_{miss}^{unused} #minus #it{p}_{miss}^{kin. fit} (GeV/#it{c})",                 binning = (600, -6, 6),     **kwargs)
    overlayCases(dataSample, "UnusedDeltaPOverP", axisTitles = "(#it{p}_{miss}^{unused} #minus #it{p}_{miss}^{kin. fit}) / #it{p}_{miss}^{kin. fit}", binning = (500, -2, 2),     **kwargs)
    overlayCases(dataSample, "UnusedDeltaTheta",  axisTitles = "#it{#theta}_{miss}^{unused} #minus #it{#theta}_{miss}^{kin. fit} (deg)",              binning = (200, -100, 100), **kwargs)
    overlayCases(dataSample, "UnusedDeltaPhi",    axisTitles = "#it{#phi}_{miss}^{unused} #minus #it{#phi}_{miss}^{kin. fit} (deg)",                  binning = (360, -180, 180), **kwargs)
    # unusedTrackData = dataSample.Filter("(NmbUnusedTracks == 1)")  # make sure unused track info exists; NOTE! this assumes that there is maximum 1 unused track
    # plot2D(unusedTrackData, xVariable = ("UnusedP_",     "UnusedP[0]"),     yVariable = "MissingProtonP",     axisTitles = "#it{p}_{miss}^{unused} (GeV/#it{c});#it{p}_{miss}^{kin. fit} (GeV/#it{c})", binning = (400, 0, 9, 400, 0, 9),           **kwargs)
    # plot2D(unusedTrackData, xVariable = ("UnusedTheta_", "UnusedTheta[0]"), yVariable = "MissingProtonTheta", axisTitles = "#it{#theta}_{miss}^{unused} (deg);#it{#theta}_{miss}^{kin. fit} (deg)",     binning = (360, 0, 180, 360, 0, 180),       **kwargs)
    # plot2D(unusedTrackData, xVariable = ("UnusedPhi_",   "UnusedPhi[0]"),   yVariable = "MissingProtonPhi",   axisTitles = "#it{#phi}_{miss}^{unused} (deg);#it{#phi}_{miss}^{kin. fit} (deg)",         binning = (360, -180, 180, 360, -180, 180), **kwargs)


if __name__ == "__main__":
  #TODO add command-line interface
  printGitInfo()
  ROOT.gROOT.SetBatch(True)
  #TODO cannot change multithreading after data frame was instantiated
  # ROOT.EnableImplicitMT(20)  # activate implicit multi-threading for RDataFrame; disable using ROOT.DisableImplicitMT()
  setupPlotStyle()

  dataPeriods = [
    "2017_01-ver03",
    "2018_01-ver02",
    "2018_08-ver02",
    "2019_11-ver01",
  ]
  treeName = "pippippimpimpmiss"
  pdfBaseDirName = "./plots"
  maxNmbTopologies = 10

  # open input files with histograms
  dataSamplesToOverlay: dict[str, dict[str, Any]] = {}
  for index, dataPeriod in enumerate(dataPeriods):
    inputFileName = f"./data/MCbggen/{dataPeriod}/{treeName}.MCbggen_{dataPeriod}.root"
    print(f"Reading generated MC histograms from file {inputFileName}")
    dataSamplesToOverlay[dataPeriod] = {
      "TFile"        : ROOT.TFile.Open(inputFileName, "READ"),
      "normToThis"   : True if index == 0 else False,
      # define plot style
      "SetLineColor" : getCbFriendlyRootColor(index, skipBlack = True),
      "SetLineWidth" : 2,
    }
  # plot generated MC truth for signal process
  histInfos = (
    {
      "histNameSuffix" : "",
      "histTitle"      : "",
    },
    {
      "histNameSuffix" : "_BeamEnergyRange1",
      "histTitle"      : "5.5 < #it{E}_{beam} < 11.0 GeV",
    },
    {
      "histNameSuffix" : "_BeamEnergyRange2",
      "histTitle"      : "7.5 < #it{E}_{beam} < 9.0 GeV",
    },
  )
  for histInfo in histInfos:
    kwargs = {
      "pdfFileNameSuffix" : "_SigMcTruth",
      "pdfDirName"        : makeDirPath(f"{pdfBaseDirName}/MCbggen"),
      "histTitle"         : histInfo["histTitle"],
    }
    overlayDataSamples1D(dataSamplesToOverlay, histName = f"SignalTruthBeamEnergy{histInfo  ['histNameSuffix']}", axisTitles = "#it{E}_{beam}^{truth} (GeV)",          **kwargs)
    overlayDataSamples1D(dataSamplesToOverlay, histName = f"SignalTruthProtonP{histInfo     ['histNameSuffix']}", axisTitles = "#it{p}_{#it{p}}^{truth} (GeV/#it{c})", **kwargs)
    overlayDataSamples1D(dataSamplesToOverlay, histName = f"SignalTruthProtonTheta{histInfo ['histNameSuffix']}", axisTitles = "#it{#theta}_{#it{p}}^{truth} (deg)",   **kwargs)
    overlayDataSamples1D(dataSamplesToOverlay, histName = f"SignalTruthProtonPhi{histInfo   ['histNameSuffix']}", axisTitles = "#it{#phi}_{#it{p}}^{truth} (deg)",     **kwargs)
    overlayDataSamples1D(dataSamplesToOverlay, histName = f"SignalTruthFourPionMass{histInfo['histNameSuffix']}", axisTitles = "#it{m}_{#it{#pi}^{#plus}#it{#pi}^{#minus}#it{#pi}^{#plus}#it{#pi}^{#minus}}^{truth} (GeV/#it{c}^{2})", **kwargs)

  # open input files with trees
  inputData: dict[str, dict[str, ROOT.RDataFrame]] = defaultdict(dict)  # dict[<data period>][<data type>]
  for dataPeriod in dataPeriods:
    inputFileNames = {
      "MCbggen" : f"./data/MCbggen/{dataPeriod}/{treeName}_flatTree.MCbggen_{dataPeriod}.root",
      "RD"      : f"./data/RD/{dataPeriod}/{treeName}_flatTree.RD_{dataPeriod}_*.root",
    }
    print(f"Reading tree '{treeName}' from files {inputFileNames}")
    for dataType, inputFileName in inputFileNames.items():
      inputData[dataPeriod][dataType] = ROOT.RDataFrame(treeName, inputFileName) \
                                            .Define("TrackFound", UNUSED_TRACK_FOUND_CONDITION) \
                                            .Filter("(-0.25 < MissingMassSquared_Measured) and (MissingMassSquared_Measured < 3.75)")  # limit data to fit range

  # overlay resolutions of kinematic variables from MC truth
  diffVariableInfos: tuple[tuple[str, str, tuple[int, float, float], tuple[float, float], str], ...] = (
    ("MissingProtonP",     "TruthDeltaP",     (400,   -4,   +4), (0,  0.7), "#it{p}_{miss}^{truth} #minus #it{p}_{miss}^{kin. fit} (GeV/#it{c})"),
    ("MissingProtonTheta", "TruthDeltaTheta", (200,  -60,  +60), (0, 20  ), "#it{#theta}_{miss}^{truth} #minus #it{#theta}_{miss}^{kin. fit} (deg)"),
    ("MissingProtonPhi",   "TruthDeltaPhi",   (360, -180, +180), (0, 70  ), "#it{#phi}_{miss}^{truth} #minus #it{#phi}_{miss}^{kin. fit} (deg)"),
  )
  for diffVariableInfo in diffVariableInfos:
    args                  = diffVariableInfo[:3]
    resPlotRange          = diffVariableInfo[3]
    diffVariableAxisTitle = diffVariableInfo[4]
    dataToOverlay: list[tuple[str, ROOT.RDataFrame]] = []
    for dataPeriod, data in inputData.items():
      dataToOverlay.append((dataPeriod, data["MCbggen"]))
    kwargs = {
      "pdfDirName"            : makeDirPath("./plots/MCbggen"),
      "additionalFilter"      : '(NmbUnusedShowers == 0) && (NmbTruthTracks == 1) && (ThrownTopology.GetString() == "2#pi^{#plus}2#pi^{#minus}p")',
      "diffVariableAxisTitle" : diffVariableAxisTitle,
    }
    overlayResolutions(dataToOverlay, *args, resBinningVariable = "BeamEnergy",         resBinning = ( 90,    2.9, 11.9), resPlotRange = resPlotRange, **kwargs)
    overlayResolutions(dataToOverlay, *args, resBinningVariable = "MissingProtonP",     resBinning = (100,    0,    5  ), resPlotRange = resPlotRange, **kwargs)
    overlayResolutions(dataToOverlay, *args, resBinningVariable = "MissingProtonTheta", resBinning = ( 56,    0,   70  ), resPlotRange = resPlotRange, **kwargs)
    overlayResolutions(dataToOverlay, *args, resBinningVariable = "MissingProtonPhi",   resBinning = ( 72, -180, +180  ), resPlotRange = resPlotRange, **kwargs)

  # overlay resolutions of kinematic variables from unused tracks
  diffVariableInfos = (
    ("MissingProtonP",     "UnusedDeltaP",     (400,   -4,   +4), (0,   1), "#it{p}_{miss}^{unused} #minus #it{p}_{miss}^{kin. fit} (GeV/#it{c})"),
    ("MissingProtonTheta", "UnusedDeltaTheta", (200,  -60,  +60), (0,  30), "#it{#theta}_{miss}^{unused} #minus #it{#theta}_{miss}^{kin. fit} (deg)"),
    ("MissingProtonPhi",   "UnusedDeltaPhi",   (360, -180, +180), (0, 110), "#it{#phi}_{miss}^{unused} #minus #it{#phi}_{miss}^{kin. fit} (deg)"),
  )
  for diffVariableInfo in diffVariableInfos:
    args                  = diffVariableInfo[:3]
    resPlotRange          = diffVariableInfo[3]
    diffVariableAxisTitle = diffVariableInfo[4]
    for dataType in ("MCbggen", "RD"):
      dataToOverlay = []
      for dataPeriod, data in inputData.items():
        dataToOverlay.append((dataPeriod, data[dataType]))
      kwargs = {
        "pdfDirName"            : makeDirPath(f"./plots/{dataType}"),
        "additionalFilter"      : "(NmbUnusedShowers == 0) && (NmbUnusedTracks == 1)",
        "diffVariableAxisTitle" : diffVariableAxisTitle,
        "pdfFileNameSuffix"     : "_unused",
      }
      overlayResolutions(dataToOverlay, *args, resBinningVariable = "BeamEnergy",         resBinning = ( 90,    2.9, 11.9), resPlotRange = resPlotRange, **kwargs)
      overlayResolutions(dataToOverlay, *args, resBinningVariable = "MissingProtonP",     resBinning = (100,    0,    5  ), resPlotRange = resPlotRange, **kwargs)
      overlayResolutions(dataToOverlay, *args, resBinningVariable = "MissingProtonTheta", resBinning = ( 56,    0,   70  ), resPlotRange = resPlotRange, **kwargs)
      overlayResolutions(dataToOverlay, *args, resBinningVariable = "MissingProtonPhi",   resBinning = ( 72, -180, +180  ), resPlotRange = resPlotRange, **kwargs)

  # overlay all periods for bggen MC and real data
  dataSamplesToOverlay = {}
  if len(inputData) > 1:
    for dataType in ("MCbggen", "RD"):
      dataSamplesToOverlay = {}
      for index, dataPeriod in enumerate(inputData.keys()):
        dataSamplesToOverlay[dataPeriod] = {
          "RDataFrame"   : inputData[dataPeriod][dataType],
          "normToThis"   : True if index == 0 else False,
          # define plot style
          "SetLineColor" : getCbFriendlyRootColor(index, skipBlack = True),
          "SetLineWidth" : 2,
        }
      makeKinematicPlotsOverlays(dataSamplesToOverlay, pdfDirName = makeDirPath(f"{pdfBaseDirName}/{dataType}"))
  # overlay bggen MC and real data for each period
  for dataPeriod in inputData.keys():
    dataSamplesToOverlay = {
      "bggen MC (scaled)" : {
        "RDataFrame"   : inputData[dataPeriod]["MCbggen"],
        # define plot style
        "SetLineColor" : ROOT.kGray,
        "SetFillColor" : ROOT.kGray,
      },
      "Real Data" : {
        "RDataFrame" : inputData[dataPeriod]["RD"],
        "normToThis" : True,
      },
    }
    makeKinematicPlotsOverlays(dataSamplesToOverlay, pdfDirName = makeDirPath(f"{pdfBaseDirName}/{dataPeriod}"))

  # make Monte Carlo plots for each period
  for dataPeriod in inputData.keys():
    makeKinematicPlotsMc(inputData[dataPeriod]["MCbggen"], isMcBggen = True, pdfDirName = makeDirPath(f"{pdfBaseDirName}/MCbggen/{dataPeriod}"))

  # make general plots for each data type and period
  for dataType in ("MCbggen", "RD"):
    for dataPeriod in inputData.keys():
      makeKinematicPlotsData(inputData[dataPeriod][dataType], pdfDirName = makeDirPath(f"{pdfBaseDirName}/{dataType}/{dataPeriod}"))
