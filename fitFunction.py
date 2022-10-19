import ROOT
ROOT.gStyle.SetOptFit(True)


# see https://root.cern/doc/master/Minuit2Minimizer_8h_source.html
# * minimizer status:
#   * status = 1 : Covariance was made pos defined
#   * status = 2 : Hesse is invalid
#   * status = 3 : Edm is above max
#   * status = 4 : Reached call limit
#   * status = 5 : Any other failure
# * status of the covariance matrix
#   * status = -1 : not available (inversion failed or Hesse failed)
#   * status =  0 : available but not positive defined
#   * status =  1 : covariance only approximate
#   * status =  2 : full matrix but forced pos def
#   * status =  3 : full accurate matrix
COV_MATRIX_STATUS_CODE = {
  -1 : ("failed",           "not available (inversion failed or Hesse failed)"),
  0  : ("not pos. def.",    "available but not positive definite"),
  1  : ("approximate",      "covariance only approximate"),
  2  : ("forced pos. def.", "full matrix but forced pos definite"),
  3  : ("accurate",         "full accurate matrix")
}
# * MINOS status code of last Minos run; minimizerStatus += 10 * minosStatus
#   * `status & 1 > 0`  : invalid lower error
#   * `status & 2 > 0`  : invalid upper error
#   * `status & 4 > 0`  : invalid because maximum number of function calls exceeded
#   * `status & 8 > 0`  : a new minimum has been found
#   * `status & 16 > 0` : error is truncated because parameter is at lower/upper limit
# * HESSE status code; minimizerStatus += 100 * hesseStatus
#   * status = 1 : hesse failed
#   * status = 2 : matrix inversion failed
#   * status = 3 : matrix is not pos defined


# define the fit function and the individual components
# see https://root-forum.cern.ch/t/syntax-of-a-free-function-or-c-functor-for-tgraph-fitting/22292/3
# and https://root.cern/manual/python/#just-in-time-compilation-of-small-strings
ROOT.gROOT.LoadMacro("./doubleGaussianPol2.C++")
# define global variables because Python callables still need to exits in scope where where fit functions are drawn
doubleGaussianPol2 = ROOT.doubleGaussianPol2()
signal = ROOT.doubleGaussianPol2(ROOT.doubleGaussianPol2.signal)
gaussian1 = ROOT.doubleGaussianPol2(ROOT.doubleGaussianPol2.gaussian1)
gaussian2 = ROOT.doubleGaussianPol2(ROOT.doubleGaussianPol2.gaussian2)
background = ROOT.doubleGaussianPol2(ROOT.doubleGaussianPol2.background)


def fixZeroPar(fitFunc, parName, errToValThr = 3):
  parVal = fitFunc.GetParameter(parName)
  parErr = fitFunc.GetParError(fitFunc.GetParNumber(parName))
  if ((parVal == 0) or (parErr / parVal < errToValThr)):
    print(f"Fixing near-zero parameter '{parName}' = {parVal} +- {parErr} to zero.")
    fitFunc.FixParameter(fitFunc.GetParNumber(parName), 0)


# fit function to measured distribution
def fitDistribution(hist, particle, fitRange = None, forceCommonGaussianMean = False, fitMissingMassSquared = True):
  print(f"Fitting histogram '{hist.GetName()}', '{hist.GetTitle()}'.")
  # construct fit function and set start parameters
  if fitRange is None:
    fitRange = (hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
  doubleGaussianPol2._forceCommonGaussianMean = forceCommonGaussianMean
  fitFunc = ROOT.TF1("doubleGaussianPol2", doubleGaussianPol2, fitRange[0], fitRange[1], 8 if forceCommonGaussianMean else 9)
  commonParNames = ("p_{0}", "p_{1}", "p_{2}", "A", "r", "#sigma_{1}", "#sigma_{2}")
  fitParameters = ((*commonParNames, "#mu") if forceCommonGaussianMean else (*commonParNames, "#mu_{1}", "#mu_{2}"))
  fitFunc.SetParNames(*fitParameters)

  # ROOT.Math.MinimizerOptions().SetDefaultTolerance(0.001)
  # ROOT.Math.MinimizerOptions().SetMaxFunctionCalls(100000)
  # ROOT.Math.MinimizerOptions().SetDefaultPrintLevel(3)
  ROOT.TVirtualFitter.SetMaxIterations(100000)

  # first fitting stage: use single Gaussian on top of pol0
  if fitMissingMassSquared:
    if particle == "Proton":
      meanStartVal = 0.9383**2  # (proton mass)^2 [GeV^2]
      widthStartVal = 0.5       # [GeV^2]
    elif particle == "Pi-" or particle == "Pi+":
      meanStartVal = 0.1396**2  # (pion mass)^2 [GeV^2]
      widthStartVal = 0.2       # [GeV^2]
    else:
      raise ValueError(f"code cannot handle particles of type '{particle}'")
  else:  # missing mass distribution
    if particle == "Proton":
      meanStartVal = 0.9383  # proton mass [GeV]
      widthStartVal = 0.25    # [GeV]
    elif particle == "Pi-" or particle == "Pi+":
      meanStartVal = 0.1396  # pion mass [GeV]
      widthStartVal = 0.2    # [GeV]
    else:
      raise ValueError(f"code cannot handle particles of type '{particle}'")
  fitFunc.SetParameter("A", hist.Integral(hist.FindBin(fitRange[0]), hist.FindBin(fitRange[1]), "WIDTH"))
  fitFunc.SetParLimits(fitFunc.GetParNumber("A"), 0, 2 * fitFunc.GetParameter("A"))  # ensure positive parameter value
  fitFunc.FixParameter(fitFunc.GetParNumber("r"), 0)
  if forceCommonGaussianMean:
    fitFunc.SetParameter("#mu", meanStartVal)
  else:
    fitFunc.SetParameter("#mu_{1}", meanStartVal)
    fitFunc.FixParameter(fitFunc.GetParNumber("#mu_{2}"), meanStartVal)
  fitFunc.SetParameter("#sigma_{1}", widthStartVal)
  fitFunc.SetParLimits(fitFunc.GetParNumber("#sigma_{1}"), 0, 10)  # ensure positive parameter value
  fitFunc.FixParameter(fitFunc.GetParNumber("#sigma_{2}"), widthStartVal)
  fitFunc.SetParameter("p_{0}", 0)
  fitFunc.FixParameter(fitFunc.GetParNumber("p_{1}"), 0)
  fitFunc.FixParameter(fitFunc.GetParNumber("p_{2}"), 0)
  # hist.Fit(fitFunc, "WLRQN")
  hist.Fit(fitFunc, "RQN")

  # second fitting stage: use double Gaussian on top of pol0
  fitFunc.ReleaseParameter(fitFunc.GetParNumber("r"))
  fitFunc.SetParameter("r", 0.75)  # nudge it, otherwise it does not move from zero
  if not forceCommonGaussianMean:
    fitFunc.ReleaseParameter(fitFunc.GetParNumber("#mu_{2}"))
  fitFunc.ReleaseParameter(fitFunc.GetParNumber("#sigma_{2}"))
  fitFunc.SetParLimits(fitFunc.GetParNumber("#sigma_{2}"), 0, 10)  # ensure positive parameter value
  # hist.Fit(fitFunc, "WLRQN")
  hist.Fit(fitFunc, "RQN")

  # third fitting stage: use double Gaussian on top of pol2
  # fitFunc.ReleaseParameter(fitFunc.GetParNumber("A"))
  # fitFunc.ReleaseParameter(fitFunc.GetParNumber("#sigma_{1}"))
  # fitFunc.ReleaseParameter(fitFunc.GetParNumber("#sigma_{2}"))
  fitFunc.ReleaseParameter(fitFunc.GetParNumber("p_{1}"))
  fitFunc.ReleaseParameter(fitFunc.GetParNumber("p_{2}"))
  # fitResult = hist.Fit(fitFunc, "WLRIMSN")
  # fitResult = hist.Fit(fitFunc, "REIMSN")
  fitResult = hist.Fit(fitFunc, "REMSN")

  # fit recovery procedure
  maxNmbRefitAttempts = 5
  nmbRefitAttempts = 0
  while ((not fitResult.IsValid() or COV_MATRIX_STATUS_CODE[fitResult.CovMatrixStatus()][0] != "accurate")
          and nmbRefitAttempts < maxNmbRefitAttempts):
    nmbRefitAttempts += 1
    print(f"Fit did not converge. Performing refit attempt #{nmbRefitAttempts} of {maxNmbRefitAttempts}.")
    zeroParNames = ["p_{0}", "p_{1}", "p_{2}"]
    for parName in zeroParNames:
      fixZeroPar(fitFunc, parName)
    # fitResult = hist.Fit(fitFunc, "WLRIMSN")
    fitResult = hist.Fit(fitFunc, "REMSN")

  # add components of fit model to LoF of histogram
  fitFunc.SetLineColor(ROOT.kRed + 1)
  fitFunc.SetLineWidth(3)
  fitFunc.SetNpx(1000)
  hist.GetListOfFunctions().Add(fitFunc)
  signal._forceCommonGaussianMean = forceCommonGaussianMean
  sigFunc = ROOT.TF1("signal", signal, fitRange[0], fitRange[1], 9)
  sigFunc.SetLineColor(ROOT.kGreen + 2)
  gaussian1._forceCommonGaussianMean = forceCommonGaussianMean
  gauss1Func = ROOT.TF1("gaussian1", gaussian1, fitRange[0], fitRange[1], 9)
  gauss1Func.SetLineStyle(ROOT.kDashed)
  gauss1Func.SetLineColor(ROOT.kGreen + 2)
  gaussian2._forceCommonGaussianMean = forceCommonGaussianMean
  gauss2Func = ROOT.TF1("gaussian2", gaussian2, fitRange[0], fitRange[1], 9)
  gauss2Func.SetLineStyle(ROOT.kDashed)
  gauss2Func.SetLineColor(ROOT.kGreen + 2)
  background._forceCommonGaussianMean = forceCommonGaussianMean
  bgFunc = ROOT.TF1("background", background, fitRange[0], fitRange[1], 9)
  bgFunc.SetLineColor(ROOT.kBlue)
  components = (sigFunc, gauss1Func, gauss2Func, bgFunc)
  for func in components:
    func.SetLineWidth(2)
    func.SetNpx(1000)
    func.SetParNames(*fitParameters)
    func.SetParameters(fitFunc.GetParameters())
    hist.GetListOfFunctions().Add(func)

  fitResult.Print()
  print(f"    reduced chi^2 = {fitResult.Chi2() / fitResult.Ndf()}; P-value = {fitResult.Prob()}")
  return fitResult
