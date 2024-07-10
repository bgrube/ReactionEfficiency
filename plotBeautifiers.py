from __future__ import annotations

from collections.abc import (
  Generator,
  Sequence,
)
from contextlib import contextmanager
from dataclasses import dataclass, field
from enum import Enum
import functools
import numbers
from typing import (
  Optional,
  Union,
)

import ROOT


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


@contextmanager
def padOf(obj: ROOT.TObject) -> Generator[Optional[ROOT.TVirtualPad], None, None]:
  """Context manager that sets pad, which contains the given TObject, as current pad"""
  lastgPad = ROOT.gPad;  # save current gPad
  # find pad of obj
  objPad: Optional[ROOT.TVirtualPad] = None
  for canv in ROOT.gROOT.GetListOfCanvases():
    # search in canvas
    if canv.GetListOfPrimitives().FindObject(obj):
      objPad = canv
      break
    # search in subpads (this probably does not work for deeper nested canvases)
    for primitive in canv.GetListOfPrimitives():
      if primitive.InheritsFrom(ROOT.TVirtualPad.Class()) and primitive.GetListOfPrimitives().FindObject(obj):
        objPad = primitive
  try:
    if objPad is not None:
      objPad.cd()
    yield objPad
  finally:
    lastgPad.cd()


@dataclass
class Lines:
  """Draws horizontal or vertical lines at given user coordinates into 1D content objects"""
  Orientation = Enum("Orientation", ("horizontal", "vertical"))

  defaultColor:         ROOT.Color_t = ROOT.kRed + 1
  defaultStyle:         ROOT.Style_t = ROOT.kDashed
  orientation:          Orientation  = Orientation.vertical
  drawContentOverLines: bool         = False
  _lineDefs:            list[tuple[float, ROOT.Color_t, ROOT.Style_t]] = field(default_factory = list)

  def set(
    self,
    lineDefs: Union[Sequence[float], Sequence[tuple[float, ROOT.Color_t, ROOT.Style_t]]]
  ) -> Lines:
    """Returns copy of object with a line defined for each given user coordinate, optionally with individual style definitions"""
    assert len(lineDefs) > 0, f"set() must be called with non-empty sequence"
    return Lines(
      defaultColor         = self.defaultColor,
      defaultStyle         = self.defaultStyle,
      orientation          = self.orientation,
      drawContentOverLines = self.drawContentOverLines,
      _lineDefs            = ([(pos, self.defaultColor, self.defaultStyle) for pos in lineDefs] if isinstance(lineDefs[0], numbers.Real)
                              else [lineDef for lineDef in lineDefs]),
    )

  def draw(
    self,
    obj: ROOT.TObject,
  ) -> None:
    """Draws the lines into the pad that contains the given TObject, e.g. TH1 or TGraph"""
    with padOf(obj) as pad:
      if pad is not None:
        line = ROOT.TLine()
        for pos, color, style in self._lineDefs:
          line.SetLineColor(color)
          line.SetLineStyle(style)
          pad.Update()
          if self.orientation == self.Orientation.vertical:
            line.DrawLine(pos, pad.GetUymin(), pos, pad.GetUymax())
          elif self.orientation == self.Orientation.horizontal:
            line.DrawLine(pad.GetUxmin(), pos, pad.GetUxmax(), pos)
        contentTypes = (ROOT.TH1, ROOT.THStack, ROOT.TGraph, ROOT.TMultiGraph)
        if (self.drawContentOverLines and any((isinstance(obj, contentType) for contentType in contentTypes))):
          # paint object over the lines
          drawOptions = obj.GetDrawOption().upper()
          if isinstance(obj, ROOT.TH1) or isinstance(obj, ROOT.THStack):
            drawOptions = f"{drawOptions} SAME"
          if isinstance(obj, ROOT.TGraph) or isinstance(obj, ROOT.TMultiGraph):
            drawOptions = drawOptions.replace("A", "")  # !Note! graphs do not have SAME option
          obj.Draw(drawOptions)
        pad.RedrawAxis()


#TODO define zero-line beautifier
