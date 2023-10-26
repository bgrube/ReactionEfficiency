import ctypes
from enum import Enum
import functools
import numpy as np
from typing import (
  Any,
  Dict,
  List,
  Optional,
  Sequence,
  Set,
  Tuple,
)

from uncertainties import UFloat, ufloat

import ROOT


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


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

def getCbFriendlyRootColor(
  index:     int,
  skipBlack: bool = False,  # if set black color is not used
) -> int:
  """Returns ROOT color index for given index in colorblind-friendly palette"""
  return getRootColor(COLORS_CB_FRIENDLY[index + (1 if skipBlack else 0)])

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
  #TODO remove dependency from external file or add file to repo
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


def calcRatioOfGraphs(graphs: Sequence[ROOT.TGraphErrors]) -> ROOT.TGraphErrors:
  """Creates graph with ratio graphs[0] / graphs[1] for points with identical x positions"""
  assert len(graphs) == 2, f"Need exactly 2 graphs to calculate ratio but got {graphs}"
  xVals: Tuple[Tuple[float, ...], Tuple[float, ...]] = (tuple(graphs[0].GetX()),  tuple(graphs[1].GetX()) )
  xErrs: Tuple[Tuple[float, ...], Tuple[float, ...]] = (tuple(graphs[0].GetEX()), tuple(graphs[1].GetEX()))
  yVals: Tuple[Tuple[float, ...], Tuple[float, ...]] = (tuple(graphs[0].GetY()),  tuple(graphs[1].GetY()) )
  yErrs: Tuple[Tuple[float, ...], Tuple[float, ...]] = (tuple(graphs[0].GetEY()), tuple(graphs[1].GetEY()))
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
      if xVals[1][j] == xVals[0][i] and xErrs[1][j] == xErrs[0][i]:
        yVal1 = ufloat(yVals[1][j], yErrs[1][j])
        break
    if yVal1 is None:
      continue
    # calculate ratio of y values and add point to ratio graph
    ratio = yVal0 / yVal1
    ratioGraph.SetPoint     (countPoints, xVal0.nominal_value, ratio.nominal_value)
    ratioGraph.SetPointError(countPoints, xVal0.std_dev,       ratio.std_dev)
    countPoints += 1
  return ratioGraph


def getGraph1DFromValues(
  graphValues:     Sequence[Tuple[UFloat, UFloat]],
  shiftByFraction: float = 0,
) -> ROOT.TGraphErrors:
  """Creates ROOT.TGraphErrors from given values"""
  if not graphValues:
    print("No data to plot")
    return
  xVals = np.array([graphVal[0].nominal_value for graphVal in graphValues], dtype = "d")
  xErrs = np.array([graphVal[0].std_dev       for graphVal in graphValues], dtype = "d")
  yVals = np.array([graphVal[1].nominal_value for graphVal in graphValues], dtype = "d")
  yErrs = np.array([graphVal[1].std_dev       for graphVal in graphValues], dtype = "d")
  # shift x values by fraction of total x range
  if shiftByFraction != 0:
    xRange = max(xVals) - min(xVals)
    shift  = xRange * shiftByFraction
    xVals = xVals + shift
  # report weighted average
  meanVal = np.average(yVals, weights = [1 / (yErr**2) for yErr in yErrs])
  print(f"    weighted mean = {meanVal}")
  return ROOT.TGraphErrors(len(xVals), xVals, yVals, xErrs, yErrs)


def getGraph2DFromValues(graphValues: Sequence[Tuple[UFloat, UFloat, UFloat]]) -> Optional[ROOT.TGraph2DErrors]:
  """Creates ROOT.TGraph2DErrors from given values"""
  if not graphValues:
    print("No data to plot")
    return None
  xVals = np.array([graphVal[0].nominal_value for graphVal in graphValues], dtype = "d")
  xErrs = np.array([graphVal[0].std_dev       for graphVal in graphValues], dtype = "d")
  yVals = np.array([graphVal[1].nominal_value for graphVal in graphValues], dtype = "d")
  yErrs = np.array([graphVal[1].std_dev       for graphVal in graphValues], dtype = "d")
  zVals = np.array([graphVal[2].nominal_value for graphVal in graphValues], dtype = "d")
  zErrs = np.array([graphVal[2].std_dev       for graphVal in graphValues], dtype = "d")
  # report weighted average
  meanVal = np.average(zVals, weights = [1 / (zErr**2) for zErr in zErrs])
  print(f"    weighted mean = {meanVal}")
  return ROOT.TGraph2DErrors(len(xVals), xVals, yVals, zVals, xErrs, yErrs, zErrs)
