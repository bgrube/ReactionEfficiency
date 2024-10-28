from __future__ import annotations

from collections.abc import Sequence
import ctypes
from enum import Enum
import functools
import numpy as np
import os
import subprocess
from typing import Any

from uncertainties import UFloat, ufloat

import ROOT


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


# black + 7 colorblind-friendly colors rom M. Okabe and K. Ito, "How to make figures and presentations that are friendly to color blind people," University of Tokyo, 2002.
# see also Bang Wong, https://www.nature.com/articles/nmeth.1618.pdf
#     https://davidmathlogic.com/colorblind
#     https://yoshke.org/blog/colorblind-friendly-diagrams
COLORS_CB_FRIENDLY: tuple[str, ...] = (
  "#000000",  # black
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#009E73",  # bluish green
  "#CC79A7",  # reddish purple
  "#56B4E9",  # sky blue
  "#E69F00",  # orange
  "#F0E442",  # yellow
)


def printGitInfo() -> None:
  """Prints directory of this file and git hash in this directory"""
  repoDir = os.path.dirname(os.path.abspath(__file__))
  gitInfo = subprocess.check_output(["git", "describe", "--always"], cwd = repoDir).strip().decode()
  print(f"Running code in '{repoDir}', git version '{gitInfo}'")


def makeDirPath(dirPath: str) -> str:
  """Create path to directory and return directory path as given"""
  try:
    os.makedirs(dirPath, exist_ok = False)
  except FileExistsError:
    pass  # directory already exists; do nothing
  except Exception:
    raise  # something went wrong
  else:
    print(f"Created directory '{dirPath}'")
  return dirPath


#TODO move into separate module
def getRootColor(hexColor: str) -> int:
  """Returns ROOT color index for given hex string in form #RRGGBB; if color does not exist yet in ROOT it is created"""
  ROOT.TColor.SetColorThreshold(0)  # ensure GetColor() returns exact color
  return ROOT.TColor.GetColor(hexColor)

def getCbFriendlyRootColor(
  index:     int,
  skipBlack: bool = False,  # if set black color is not used
) -> int:
  """Returns ROOT color index for given index in colorblind-friendly palette"""
  return getRootColor(COLORS_CB_FRIENDLY[index + (1 if skipBlack else 0)])

# 11 filled marker styles; the float is a relative scaling factor to obtain equal apparent sizes
MARKERS_FILLED: tuple[tuple[int, float], ...] = (
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
MARKERS_OPEN: tuple[tuple[int, float], ...] = (
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

#TODO take TObject and check wether it is TAttLine etc.
#TODO add cycle option; provide more styles by combining colors and markers
def setCbFriendlyStyle(
  graphOrHist:   ROOT.TGraph | ROOT.TH1,
  styleIndex:    int,  # index that switches between styles
  skipBlack:     bool  = True,  # if set black color is not used
  setMarker:     bool  = True,
  markerSize:    float = 1.5,
  filledMarkers: bool  = True,
) -> None:
  """Sets line color and marker style, color, and size of a TGraph or TH1 according to a style index"""
  nmbStyles = min(len(COLORS_CB_FRIENDLY) - (1 if skipBlack else 0), len(MARKERS_FILLED), len(MARKERS_OPEN))
  assert styleIndex < nmbStyles, f"The style index {styleIndex} goes beyond the maximum of {nmbStyles} styles that are implemented"
  color = getCbFriendlyRootColor(styleIndex, skipBlack)
  graphOrHist.SetLineColor(color)
  if setMarker:
    graphOrHist.SetMarkerColor(color)
    graphOrHist.SetMarkerStyle(MARKERS_FILLED[styleIndex][0] if filledMarkers else MARKERS_OPEN[styleIndex][0])
    graphOrHist.SetMarkerSize(markerSize * MARKERS_FILLED[styleIndex][1] if filledMarkers else MARKERS_OPEN[styleIndex][1])


def setupPlotStyle(rootLogonPath: str = "./rootlogon.C") -> None:
  """Defines ROOT plot style"""
  ROOT.gROOT.LoadMacro(rootLogonPath)
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


def getRangeOfGraph(graph: ROOT.TGraph) -> tuple[float, float, float, float]:
  xMin = ctypes.c_double()
  xMax = ctypes.c_double()
  yMin = ctypes.c_double()
  yMax = ctypes.c_double()
  graph.ComputeRange(xMin, yMin, xMax, yMax)
  return (xMin.value, yMin.value, xMax.value, yMax.value)


#TODO replace by Lines beautifier
def drawZeroLine(
  obj: ROOT.TObject,
  style = ROOT.kDashed,
  color = ROOT.kBlack
) -> None:
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


def redrawFrame(pad: ROOT.TVirtualPad) -> None:
  """redraws histogram frame to fix overprinting by histogram content"""
  # unfortunately, filled TH1 or TH2 drawn with COLZ are drawn over
  # the histogram frame so that they mask half of the line width of
  # the frame
  # even worse, this is not considered a bug: https://root-forum.cern.ch/t/th2-colz-obscures-the-top-frame/16699
  # the "official" workaround is to redraw the frame boc by hand; sigh!
  pad.RedrawAxis()
  pad.Update()
  frame = ROOT.TBox()
  frame.SetLineColor(ROOT.gStyle.GetFrameLineColor())
  frame.SetLineStyle(ROOT.gStyle.GetFrameLineStyle())
  frame.SetLineWidth(ROOT.gStyle.GetFrameLineWidth())
  frame.SetFillStyle(0)
  xMin = pad.GetUxmin()
  xMax = pad.GetUxmax()
  yMin = pad.GetUymin()
  yMax = pad.GetUymax()
  if (pad.GetLogx() == 1):
    xMin = pow(10, xMin)
    xMax = pow(10, xMax)
  if (pad.GetLogy() == 1):
    yMin = pow(10, yMin)
    yMax = pow(10, yMax)
  frame.DrawBox(xMin, yMin, xMax, yMax)


def callMemberFunctionsWithArgs(
  instance:          Any,             # instance for which member functions will be called
  functionsWithArgs: dict[str, Any],  # member-function names with arguments
) -> None:
  """Calls member functions of given object with given arguments"""
  for functionName, argument in functionsWithArgs.items():
    function = getattr(instance, functionName, None)
    if function is None:
      continue
    # print(f"Calling member function '{functionName}({argument})' of {instance}")
    function(argument)


def calcRatioOfGraphs1D(
  graphs:     Sequence[ROOT.TGraphErrors],
  ratioRange: tuple[float | None, float | None] = (None, None),  # is set, ratios outside this range are not filled into graph
) -> ROOT.TGraphErrors:
  """Creates 1D graph with ratio graphs[0] / graphs[1] for points with identical x positions"""
  assert len(graphs) == 2, f"Need exactly 2 graphs to calculate ratio but got {graphs}"
  xVals: tuple[tuple[float, ...], tuple[float, ...]] = (tuple(graphs[0].GetX()),  tuple(graphs[1].GetX()) )
  xErrs: tuple[tuple[float, ...], tuple[float, ...]] = (tuple(graphs[0].GetEX()), tuple(graphs[1].GetEX()))
  yVals: tuple[tuple[float, ...], tuple[float, ...]] = (tuple(graphs[0].GetY()),  tuple(graphs[1].GetY()) )
  yErrs: tuple[tuple[float, ...], tuple[float, ...]] = (tuple(graphs[0].GetEY()), tuple(graphs[1].GetEY()))
  ratioGraph = ROOT.TGraphErrors()
  ratioGraph.SetName(f"{graphs[0].GetName()}_{graphs[1].GetName()}")
  countPoints = 0  # counts points that match in both graphs
  # loop over data points in first graph
  for i in range(graphs[0].GetN()):
    xVal0 = ufloat(xVals[0][i], xErrs[0][i])
    yVal0 = ufloat(yVals[0][i], yErrs[0][i])
    # find matching data point with same x position in second graph
    yVal1 = None
    for j in range(graphs[1].GetN()):
      if xVals[1][j] == xVals[0][i] and xErrs[1][j] == xErrs[0][i]:  #TODO check whether there are rounding issues
        yVal1 = ufloat(yVals[1][j], yErrs[1][j])
        break
    if yVal1 is None:
      continue
    # calculate ratio of y values and add point to ratio graph
    ratio = yVal0 / yVal1
    if ratioRange[0] is not None and ratio.nominal_value < ratioRange[0]:
      continue
    if ratioRange[1] is not None and ratio.nominal_value > ratioRange[1]:
      continue
    ratioGraph.SetPoint     (countPoints, xVal0.nominal_value, ratio.nominal_value)
    ratioGraph.SetPointError(countPoints, xVal0.std_dev,       ratio.std_dev)
    countPoints += 1
  return ratioGraph


def calcRatioOfGraphs2D(
  graphs:     Sequence[ROOT.TGraph2DErrors],
  ratioRange: tuple[float | None, float | None] = (None, None),  # is set, ratios outside this range are not filled into graph
) -> ROOT.TGraph2DErrors:
  """Creates 2D graph with ratio graphs[0] / graphs[1] for points with identical (x, y) positions"""
  assert len(graphs) == 2, f"Need exactly 2 graphs to calculate ratio but got {graphs}"
  xVals: tuple[tuple[float, ...], tuple[float, ...]] = (tuple(graphs[0].GetX()),  tuple(graphs[1].GetX()) )
  xErrs: tuple[tuple[float, ...], tuple[float, ...]] = (tuple(graphs[0].GetEX()), tuple(graphs[1].GetEX()))
  yVals: tuple[tuple[float, ...], tuple[float, ...]] = (tuple(graphs[0].GetY()),  tuple(graphs[1].GetY()) )
  yErrs: tuple[tuple[float, ...], tuple[float, ...]] = (tuple(graphs[0].GetEY()), tuple(graphs[1].GetEY()))
  zVals: tuple[tuple[float, ...], tuple[float, ...]] = (tuple(graphs[0].GetZ()),  tuple(graphs[1].GetZ()) )
  zErrs: tuple[tuple[float, ...], tuple[float, ...]] = (tuple(graphs[0].GetEZ()), tuple(graphs[1].GetEZ()))
  ratioGraph = ROOT.TGraph2DErrors()
  ratioGraph.SetName(f"{graphs[0].GetName()}_{graphs[1].GetName()}")
  countPoints = 0  # counts points that match in both graphs
  # loop over data points in first graph
  for i in range(graphs[0].GetN()):
    xVal0 = ufloat(xVals[0][i], xErrs[0][i])
    yVal0 = ufloat(yVals[0][i], yErrs[0][i])
    zVal0 = ufloat(zVals[0][i], zErrs[0][i])
    # find matching data point with same (x, y) position in second graph
    zVal1 = None
    for j in range(graphs[1].GetN()):
      if (xVals[1][j] == xVals[0][i] and xErrs[1][j] == xErrs[0][i]) and (yVals[1][j] == yVals[0][i] and yErrs[1][j] == yErrs[0][i]):  #TODO check whether there are rounding issues
        zVal1 = ufloat(zVals[1][j], zErrs[1][j])
        break
    if zVal1 is None:
      continue
    # calculate ratio of z values and add point to ratio graph
    zRatio = zVal0 / zVal1
    if ratioRange[0] is not None and zRatio.nominal_value < ratioRange[0]:
      continue
    if ratioRange[1] is not None and zRatio.nominal_value > ratioRange[1]:
      continue
    ratioGraph.SetPoint     (countPoints, xVal0.nominal_value, yVal0.nominal_value, zRatio.nominal_value)
    ratioGraph.SetPointError(countPoints, xVal0.std_dev,       yVal0.std_dev,       zRatio.std_dev)
    countPoints += 1
  return ratioGraph


def getGraph1DFromValues(
  graphValues:     Sequence[tuple[UFloat, UFloat]],
  shiftByFraction: float = 0,
) -> ROOT.TGraphErrors | None:
  """Creates ROOT.TGraphErrors from given values"""
  if not graphValues:
    print("No data to plot")
    return None
  xVals = np.array([graphVal[0].nominal_value for graphVal in graphValues], dtype = np.float64)
  xErrs = np.array([graphVal[0].std_dev       for graphVal in graphValues], dtype = np.float64)
  yVals = np.array([graphVal[1].nominal_value for graphVal in graphValues], dtype = np.float64)
  yErrs = np.array([graphVal[1].std_dev       for graphVal in graphValues], dtype = np.float64)
  # shift x values by fraction of total x range
  if shiftByFraction != 0:
    xRange = max(xVals) - min(xVals)
    shift  = xRange * shiftByFraction
    xVals = xVals + shift
  # report weighted average
  meanVal = np.average(yVals, weights = [1 / (yErr**2) for yErr in yErrs])
  print(f"    weighted mean = {meanVal}")
  return ROOT.TGraphErrors(len(xVals), xVals, yVals, xErrs, yErrs)


def getGraph2DFromValues(graphValues: Sequence[tuple[UFloat, UFloat, UFloat]]) -> ROOT.TGraph2DErrors | None:
  """Creates ROOT.TGraph2DErrors from given values"""
  if not graphValues:
    print("No data to plot")
    return None
  xVals = np.array([graphVal[0].nominal_value for graphVal in graphValues], dtype = np.float64)
  xErrs = np.array([graphVal[0].std_dev       for graphVal in graphValues], dtype = np.float64)
  yVals = np.array([graphVal[1].nominal_value for graphVal in graphValues], dtype = np.float64)
  yErrs = np.array([graphVal[1].std_dev       for graphVal in graphValues], dtype = np.float64)
  zVals = np.array([graphVal[2].nominal_value for graphVal in graphValues], dtype = np.float64)
  zErrs = np.array([graphVal[2].std_dev       for graphVal in graphValues], dtype = np.float64)
  # report weighted average
  meanVal = np.average(zVals, weights = [1 / (zErr**2) for zErr in zErrs])
  print(f"    weighted mean = {meanVal}")
  return ROOT.TGraph2DErrors(len(xVals), xVals, yVals, zVals, xErrs, yErrs, zErrs)


Graph2DVar = Enum("Graph2DVar", ("x", "y"))

def slice2DGraph(
  graph2D:     ROOT.TGraph2DErrors,
  steppingVar: Graph2DVar,  # each 1D graph corresponds to a slice in this variable
) -> dict[tuple[float, float], ROOT.TGraphErrors]:
  """Slices a 2D graph into 1D graphs assuming equidistant binning in stepping variable; returns dictionary with bin range of stepping variable and corresponding 1D graph"""
  # read values from 2D graph assuming equidistant binning
  values: dict[str, tuple[float, ...]] = {
    "x"    : tuple(graph2D.GetX()),
    "xErr" : tuple(graph2D.GetEX()),
    "y"    : tuple(graph2D.GetY()),
    "yErr" : tuple(graph2D.GetEY()),
    "z"    : tuple(graph2D.GetZ()),
    "zErr" : tuple(graph2D.GetEZ()),
  }
  steppingVarBinCenters:    set[float] = set(values[steppingVar.name])
  steppingVarHalfBinWidths: set[float] = set((round(value, 10) for value in values[steppingVar.name + "Err"]))
  assert len(steppingVarHalfBinWidths) == 1, f"Binning for stepping variable is not equidistant; found half bin widths {steppingVarHalfBinWidths}"
  steppingVarHalfBinWidth = next(iter(steppingVarHalfBinWidths))
  plottingVar = Graph2DVar.x if steppingVar == Graph2DVar.y else Graph2DVar.y
  # construct 1D graph for each bin of stepping variable
  graphs1D: dict[tuple[float, float], ROOT.TGraphErrors] = {}
  for steppingVarBinCenter in sorted(steppingVarBinCenters):
    xVals = np.array([value for i, value in enumerate(values[plottingVar.name])         if values[steppingVar.name][i] == steppingVarBinCenter], dtype = np.float64)
    xErrs = np.array([value for i, value in enumerate(values[plottingVar.name + "Err"]) if values[steppingVar.name][i] == steppingVarBinCenter], dtype = np.float64)
    yVals = np.array([value for i, value in enumerate(values["z"])                      if values[steppingVar.name][i] == steppingVarBinCenter], dtype = np.float64)
    yErrs = np.array([value for i, value in enumerate(values["zErr"])                   if values[steppingVar.name][i] == steppingVarBinCenter], dtype = np.float64)
    assert len(xVals) > 0, f"Could not find any values in graph '{graph2D.GetName()}' for bin center {steppingVar} == {steppingVarBinCenter}"
    steppingVarBinRange = (round(steppingVarBinCenter - steppingVarHalfBinWidth, 6), round(steppingVarBinCenter + steppingVarHalfBinWidth, 6))
    graphs1D[steppingVarBinRange] = ROOT.TGraphErrors(len(xVals), xVals, yVals, xErrs, yErrs)
  return graphs1D
