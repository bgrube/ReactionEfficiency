#!/usr/bin/env python3.11

import ROOT


if __name__ == "__main__":

  # x = ROOT.RooRealVar("x", "x", -10, 10)
  # mean = ROOT.RooRealVar("mean", "mean of gaussian", 1, -10, 10)
  # sigma = ROOT.RooRealVar("sigma", "width of gaussian", 1, 0.1, 10)
  # gauss = ROOT.RooGaussian("gauss", "gaussian PDF", x, mean, sigma)

  # w = ROOT.RooWorkspace("w")
  # gauss = w.factory("Gaussian::gauss(x[-10, 10], mean[1, -10, 10], sigma[1, 0.1, 10])")
  # x = w.var("x")
  # sigma = w.var("sigma")

  # xframe = x.frame(ROOT.RooFit.Title("Gaussian pdf"))

  # gauss.plotOn(xframe)
  # sigma.setVal(3)
  # gauss.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kRed))

  w = ROOT.RooWorkspace("w")

  # "logarithmic Gaussian"
  # Eq. (9) in NIMA 441 (2000) 401 at https://doi.org/10.1016/S0168-9002(99)00992-4
  # gaussSkewed = w.factory("Novosibirsk::gaussSkewed(x[-10, 10], peak[1, -10, 10], width[1, 0.1, 10], skew[0, -1, 1])")
  gaussSkewed = w.factory(f"Novosibirsk::gaussSkewed(x[-0.25, 3.75], peak[1.1, 0, 2], width[1.0, 0.01, 2], skew[-0.889782, -2, 2])")

  # # skew normal PDF (see https://www.wikiwand.com/en/Skew_normal_distribution)
  # # gaussSkewed = w.factory("PROD::gaussSkewed("
  # #   "Gaussian::gauss(x[-10, 10], peak[1, -10, 10], width[1, 0.1, 10]),"
  # #   "EXPR::erfc('2 * TMath::Erfc(skew * ((x - peak) / width))', x, peak, width, skew[0, -10, 10])"
  # # ")")
  # gaussSkewed = w.factory("EXPR::gaussSkewed("
  #   "'2 * TMath::Gaus(x, peak, width, true)"
  #   "* TMath::Erfc(skew * ((x - peak) / width))',"
  #   "x[-10, 10], peak[1, -10, 10], width[1, 0.1, 10], skew[0, -10, 10])")

  # # exponentially modified Gaussian PDF (see https://www.wikiwand.com/en/Exponentially_modified_Gaussian_distribution)
  # # smaller skew means higher skewness
  # # max = peak + 1 / skew
  # gaussSkewed = w.factory("EXPR::gaussSkewed("
  #   "'(skew / 2) * TMath::Exp((skew / 2) * (2 * peak + skew * width * width - 2 * x))"
  #   "* TMath::Erfc((peak + skew * width * width - x) / (sqrt(2) * width))',"
  # "x[-5, 15], peak[1, -10, 10], width[1, 0.1, 10], skew[25, 0, 25])")

  x  = w.var('x')
  erfc = w.function("gaussSkewed")
  print("xvar =", x.getVal())
  print("erfc =", erfc.getVal())
  x.setVal(2)
  print("xvar =", x.getVal())
  print("erfc =", erfc.getVal())

  xframe = x.frame(ROOT.RooFit.Title("Skewed Gaussian PDF"))
  gaussSkewed.plotOn(xframe)
  skew = w.var("skew")
  skew.setVal(-2)
  gaussSkewed.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kRed))
  c = ROOT.TCanvas()
  xframe.Draw()
  c.SaveAs("./plotRooFitPdf.pdf")

  # xvar = ROOT.RooRealVar('x', 'x', 0, -10, 10)
  # peak = ROOT.RooRealVar('peak', 'peak', 1, -10, 10)
  # width = ROOT.RooRealVar('width', 'width', 1, 0.1, 10)
  # skew = ROOT.RooRealVar('skew', 'skew', 0.1, -1, 1)
  # w.Import(xvar)
  # w.Import(peak)
  # w.Import(width)
  # w.Import(skew)
  # gaussSkewed = w.factory(
  #   "EXPR::erfc('2 * x', x)"
  # # "PROD::gaussSkewed("
  # #   "Gaussian::gauss(x[-10, 10], peak[1, -10, 10], width[1, 0.1, 10]),"
  # #   "Gaussian::gauss(x[-10, 10], peak[1, -10, 10], width[1, 0.1, 10])"
  # #   # "expr::erfc('x * x', x)"
  # # ")"
  # )
  # gaussSkewed = w.factory("EXPR::erfc('x * x + 1', x)")
  # gaussSkewed = w.factory("EXPR::erfc('TMath::Erfc(skew * ((x - peak) / width))', x, peak, width, skew)")
