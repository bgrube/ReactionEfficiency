#!/usr/bin/env python3


import argparse
import functools
import os
import sys
from typing import (
  Iterable,
  List,
  Mapping,
  Optional,
  Sequence,
  Tuple,
)

import ROOT

import plotTools


# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


def andCuts(cuts: Iterable[Optional[str]]) -> str:
  """Creates cut string where given cuts are combined by via logical and, ignoring None and empty cut strings"""
  cutsWithBraces = (f"({cut})" for cut in filter(None, cuts))  # surround each cut by braces and filter out None _and_ ""
  return " && ".join(cutsWithBraces)


#TODO convert collection of functions to class with member functions
def readWeights(
  inputFileName:     str,       # name of file from which to read weights
  inputTreeName:     str,       # name of tree from which to read weights
  weightBranchName:  str,       # name of branch from which to read weights
  comboIdName:       str,       # name of branch with unique ID
  sWeightLabel:      str,       # label used when applying weights
  sWeightFileName:   str,       # ROOT file name to which weight object will be written
  sWeightObjectName: str,       # name of weight object in ROOT file
  cut:               str = "",  # optional selection cut to apply
) -> None:
  """Reads weights from input tree and writes them out as a HS::FIT::Weights object"""
  print(f"Reading weights '{weightBranchName}' from tree '{inputTreeName}' in file '{inputFileName}'"
  f", writing them to key '{sWeightObjectName}' in file '{sWeightFileName}', and assigning label '{sWeightLabel}'"
  + ("" if not cut else f" while applying cut(s) '{cut}'"))
  currentDir = ROOT.gDirectory
  inputFile = ROOT.TFile.Open(inputFileName, "READ")
  inputTree = inputFile.Get(inputTreeName)
  currentDir.cd()
  weights = ROOT.Weights(sWeightObjectName)  # set name of the Weights object
  weights.SetFile(sWeightFileName)
  weights.SetSpecies(sWeightLabel)
  weights.SetIDName(comboIdName)
  weights.WeightBySelection(inputTree, (cut or "(1)"), weightBranchName)
  weights.SortWeights()
  weights.Save()


def defineGaussianPdf(
  fitManager:  "ROOT.FitManager",
  fitVariable: str,
  pdfName:     str,
  parDefs:     Mapping[str, str],  # maps parameter names to their definition strings
) -> None:
  """Defines Gaussian PDF"""
  fitManager.SetUp().FactoryPDF(f"Gaussian::{pdfName}({fitVariable}, mean_{pdfName}[{parDefs['mean']}], width_{pdfName}[{parDefs['width']}])")


def defineDoubleGaussianPdf(
  fitManager:  "ROOT.FitManager",
  fitVariable: str,
  pdfName:     str,
  parDefs:     Mapping[str, str],  # maps parameter names to their definition strings
  commonMean:  bool = False,
) -> None:
  """Defines two types of double-Gaussian PDFs"""
  if commonMean:
    # with same mean
    fitManager.SetUp().FactoryPDF(f"SUM::{pdfName}("
      f"r_{pdfName}[{parDefs['r']}] * Gaussian::{pdfName}_N1({fitVariable}, mean_{pdfName}[{parDefs['mean']}], width1_{pdfName}[{parDefs['width1']}]),"
                                    f"Gaussian::{pdfName}_N2({fitVariable}, mean_{pdfName},                    width2_{pdfName}[{parDefs['width2']}])"
    ")")
  else:
    # with separate means
    fitManager.SetUp().FactoryPDF(f"SUM::{pdfName}("
      f"r_{pdfName}[{parDefs['r']}] * Gaussian::{pdfName}_N1({fitVariable}, mean1_{pdfName}[{parDefs['mean1']}], width1_{pdfName}[{parDefs['width1']}]),"
                                    f"Gaussian::{pdfName}_N2({fitVariable}, mean2_{pdfName}[{parDefs['mean2']}], width2_{pdfName}[{parDefs['width2']}])"
    ")")


def defineSkewedGaussianPdf(
  fitManager:  "ROOT.FitManager",
  fitVariable: str,
  pdfName:     str,
  parDefs:     Mapping[str, str],  # maps parameter names to their definition strings
  pdfType:     str = "skewNormal",
) -> None:
  """Defines several types of skewed Gaussian PDFs"""
  # various implementations of skewed Gaussians
  if pdfType == "SkewNormal":
    # skew normal PDF (see https://en.wikipedia.org/wiki/Skew_normal_distribution)
    fitManager.SetUp().FactoryPDF(f"EXPR::{pdfName}("
      f"'2 * TMath::Gaus({fitVariable}, peak_{pdfName}, width_{pdfName}, true)"
      f"* TMath::Erfc(skew_{pdfName} * (({fitVariable} - peak_{pdfName}) / width_{pdfName}))',"
      f"{fitVariable}, peak_{pdfName}[{parDefs['peak']}], width_{pdfName}[{parDefs['width']}], skew_{pdfName}[{parDefs['skew']}]"
    ")")
  elif pdfType == "ExpMod":
    # exponentially modified Gaussian PDF (see https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution)
    # smaller skew means higher skewness
    # max = peak + 1 / skew
    fitManager.SetUp().FactoryPDF(f"EXPR::{pdfName}("
      f"'(skew_{pdfName} / 2) * TMath::Exp((skew_{pdfName} / 2) * (2 * peak_{pdfName} + skew_{pdfName} * width_{pdfName} * width_{pdfName} - 2 * {fitVariable}))"
      f"* TMath::Erfc((peak_{pdfName} + skew_{pdfName} * width_{pdfName} * width_{pdfName} - {fitVariable}) / (sqrt(2) * width_{pdfName}))',"
      f"{fitVariable}, peak_{pdfName}[{parDefs['peak']}], width_{pdfName}[{parDefs['width']}], skew_{pdfName}[{parDefs['skew']}]"
    ")")
  elif pdfType == "Log":
    # "logarithmic Gaussian"
    # Eq. (9) in NIMA 441 (2000) 401 at https://doi.org/10.1016/S0168-9002(99)00992-4
    fitManager.SetUp().FactoryPDF(f"Novosibirsk::{pdfName}({fitVariable}, peak_{pdfName}[{parDefs['peak']}], width_{pdfName}[{parDefs['width']}], skew_{pdfName}[{parDefs['skew']}])")
  else:
    #TODO try https://en.wikipedia.org/wiki/Normal-exponential-gamma_distribution
    raise ValueError(f"Unknown skewed Gaussian type '{pdfType}'")


def defineHistogramPdf(
  fitManager:           "ROOT.FitManager",
  fitVariable:          str,
  pdfName:              str,
  parDefs:              Mapping[str, str],  # maps parameter names to their definition strings
  outputDirName:        str,  # name of directory where weight files are written
  templateDataFileName: str,  # name of file from which histogram is filled
  templateDataTreeName: str,  # name of tree from which histogram is filled
  weightBranchName:     str,  # name of branch from which to read weights
  comboIdName:          str,  # name of branch with unique combo ID
  cut:                  str,  # cut that is applied when filling histogram
) -> None:
  """Defines histogram-based PDF"""
  # create RF-sideband weights for histogram PDF
  rfSWeightLabel      = "RfSideband"
  rfSWeightObjectName = f"{rfSWeightLabel}Weights{pdfName}"
  rfSWeightFileName   = f"{outputDirName}/{rfSWeightObjectName}.root"
  readWeights(
    inputFileName     = templateDataFileName,
    inputTreeName     = templateDataTreeName,
    weightBranchName  = weightBranchName,
    comboIdName       = comboIdName,
    sWeightLabel      = rfSWeightLabel,
    sWeightFileName   = rfSWeightFileName,
    sWeightObjectName = rfSWeightObjectName,
    cut               = cut,
  )
  # histogram PDF with fudge parameters that allow (small) deviations from the original shape:
  #     smear = width of Gaussian the histogram is convoluted with
  #     shift = shifts histogram in x-direction, i.e. PDF(x - shift)
  #     scale = scales histogram in x-direction around its maximum value, i.e. PDF(scale * (x - maxPos) + maxPos)
  fitManager.SetUp().FactoryPDF(
    f"RooHSEventsHistPDF::{pdfName}({fitVariable}, smear_{pdfName}[{parDefs['smear']}], shift_{pdfName}[{parDefs['shift']}], scale_{pdfName}[{parDefs['scale']}])"
    + f"WEIGHTS@{rfSWeightLabel},{rfSWeightFileName},{rfSWeightObjectName}"  # apply sWeights created above; !Note! no whitespace allowed in this string
  )
  # constrain PDF fudge parameters
  # the constraints are derived from the smear, shift, and scale parameter initial values and limits, i.e.
  #     Gaussian mean  = initial value
  #     Gaussian sigma = max / 5 for smear
  #     Gaussian sigma = range / 10 for shift and scale
  pdf = fitManager.SetUp().WS().pdf(pdfName)
  if len(parDefs['smear'].split(",")) == 3:
    fitManager.SetUp().AddGausConstraint(pdf.AlphaConstraint())  # constrain smear
  if len(parDefs['shift'].split(",")) == 3:
    fitManager.SetUp().AddGausConstraint(pdf.OffConstraint())    # constrain shift
  if len(parDefs['scale'].split(",")) == 3:
    fitManager.SetUp().AddGausConstraint(pdf.ScaleConstraint())  # constrain scale
  # load data for template histograms; ensure that calling code has already defined the kinematical binnings
  print(f"Loading data for histogram PDF '{pdfName}' from tree '{templateDataTreeName}' in file '{templateDataFileName}'"
  + f" and applying weights from '{rfSWeightObjectName}' in file '{rfSWeightFileName}'")
  fitManager.LoadSimulated(templateDataTreeName, templateDataFileName, pdfName)


def defineSigPdf(
  fitManager:           "ROOT.FitManager",
  fitVariable:          str,
  pdfType:              str,  # selects type of PDF
  fixPars:              Sequence[str] = (),  # tuple with fit-parameter names to fix
  pdfName:              str = "SigPdf",
  # arguments only needed for histogram PDF
  outputDirName:        str = "",  # name of directory where weight files are written
  templateDataFileName: str = "",  # name of file from which histogram is filled
  templateDataTreeName: str = "",  # name of tree from which histogram is filled
  weightBranchName:     str = "",  # name of branch from which to read weights
  comboIdName:          str = "",  # name of branch with unique combo ID
  cut:                  str = "",  # cut that is applied when filling histogram
) -> None:
  """Defines signal PDFs of various types"""
  print(f"Defining signal PDF '{pdfName}' of type '{pdfType}'")

  pdfTypeArgs = pdfType.split("_")
  if pdfTypeArgs[0] == "Gaussian":
    defineGaussianPdf(
      fitManager, fitVariable, pdfName,
      {
        "mean"  : f"{0.9383**2}, 0, 2",
        "width" : "0.3, 0.01, 2",
      },
    )
  elif pdfTypeArgs[0] == "DoubleGaussian":
    defineDoubleGaussianPdf(
      fitManager, fitVariable, pdfName,
      {
        "r"      : "0.5, 0, 1",  # fraction of Gaussian 1
        "mean"   : f"{0.9383**2}, 0, 2",  # common-mean case
        "mean1"  :  "1.0,         0, 2",  # separate-means case
        "mean2"  : f"{0.9383**2}, 0, 2",
        "width1" : "1.0, 0.01, 2",  # wide Gaussian
        "width2" : "0.2, 0.01, 2",  # narrow Gaussian
      },
      commonMean = len(pdfTypeArgs) == 2 and pdfTypeArgs[1] == "SameMean",
    )
  elif pdfTypeArgs[0] == "Histogram":
    defineHistogramPdf(
      fitManager, fitVariable, pdfName,
      {
        "smear" : "0" if "smear" in fixPars else "0,  0,    0.5",
        "shift" : "0" if "shift" in fixPars else "0, -0.25, 0.25",
        "scale" : "1" if "scale" in fixPars else "1,  0.5,  1.5",
      },
      outputDirName,
      templateDataFileName,
      templateDataTreeName,
      weightBranchName,
      comboIdName,
      cut = andCuts((cut, "(IsSignal == 1)")),
    )
  else:
    raise ValueError(f"Unknown signal PDF type '{pdfTypeArgs[0]}'")

  pdfWeightStartVal = 1.0
  fitManager.SetUp().LoadSpeciesPDF(pdfName, pdfWeightStartVal)


def defineBkgPdf(
  fitManager:           "ROOT.FitManager",
  fitVariable:          str,
  pdfType:              str,  # selects type of PDF
  fixPars:              Sequence[str] = (),  # tuple with fit-parameter names to fix
  pdfName:              str = "BkgPdf",
  # arguments only needed for histogram PDF
  outputDirName:        str = "",  # name of directory where weight files are written
  templateDataFileName: str = "",  # name of file from which histogram is filled
  templateDataTreeName: str = "",  # name of tree from which histogram is filled
  weightBranchName:     str = "",  # name of branch from which to read weights
  comboIdName:          str = "",  # name of branch with unique combo ID
  cut:                  str = "",  # cut that is applied when filling histogram
):
  """Defines signal PDFs of various types"""
  print(f"Defining background PDF '{pdfName}' of type '{pdfType}'")

  #TODO add polynomials
  # # 2nd-order positive-definite polynomial
  # fitManager.SetUp().FactoryPDF(f"GenericPdf::BkgPdf('@1 * @1 + (@2 + @3 * @0) * (@2 + @3 * @0)', {{{fitVariable}, p0_BkgPdf[0, -100, 100], p1_BkgPdf[0, -100, 100], p2_BkgPdf[0, -100, 100]}})")
  #
  # Chebychev polynomial
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BkgPdf({fitVariable}, {{}})")
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BkgPdf({fitVariable}, {{p0_BkgPdf[0, 0, 100]}})")
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BkgPdf({fitVariable}, {{p0_BkgPdf[0, -1, 1], p1_BkgPdf[0, -1, 1]}})")
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BkgPdf({fitVariable}, {{p0_BkgPdf[0, -1, 1], p1_BkgPdf[0, -1, 1], p2_BkgPdf[0, -1, 1]}})")
  #
  # # Bernstein polynomial
  # # see https://root.cern.ch/doc/master/classRooBernstein.html
  # fitManager.SetUp().FactoryPDF(f"Bernstein::BkgPdf({fitVariable}, {{}})")
  # fitManager.SetUp().FactoryPDF(f"Bernstein::BkgPdf({fitVariable}, {{p0_BkgPdf[0, 0, 1]}})")
  # fitManager.SetUp().FactoryPDF(f"Bernstein::BkgPdf({fitVariable}, {{p0_BkgPdf[0, 0, 1], p1_BkgPdf[0, 0, 1]}})")
  # fitManager.SetUp().FactoryPDF(f"Bernstein::BkgPdf({fitVariable}, {{p0_BkgPdf[0, 0, 1], p1_BkgPdf[0, 0, 1], p2_BkgPdf[0, 0, 1]}})")
  pdfTypeArgs = pdfType.split("_")
  if pdfTypeArgs[0] == "DoubleGaussian":
    defineDoubleGaussianPdf(
      fitManager, fitVariable, pdfName,
      {
        "r"      : "0.5, 0,    1",  # fraction of Gaussian 1
        "mean"   : "1.4, 0,    2",  # common-mean case
        "mean1"  : "1.5, 0,    2",  # separate-means case
        "mean2"  : "1.2, 0,    2",
        "width1" : "1.0, 0.01, 2",  # wide Gaussian
        "width2" : "0.4, 0.01, 2",  # narrow Gaussian
      },
      commonMean = len(pdfTypeArgs) == 2 and pdfTypeArgs[1] == "SameMean",
    )
  elif pdfTypeArgs[0] == "SkewedGaussian":
    assert len(pdfTypeArgs) == 2, f"Skewed Gaussian PDF requires subtype, i.e. 'SkewedGaussian_<subType>', given PDF type is '{pdfType}'"
    pdfSubType = pdfTypeArgs[1]
    parDefs = {}
    if pdfSubType == "SkewNormal":
      parDefs = {
        "peak"  : "0.8,  0,    2",
        "width" : "0.9,  0.01, 2",
        "skew"  : "-2.5, -5,   5",
      }
    elif pdfSubType == "ExpMod":
      parDefs = {
        "peak"  : "0.8, 0,    2",
        "width" : "0.5, 0.01, 2",
        "skew"  : "2,   0,    15",
      }
    elif pdfSubType == "Log":
      parDefs = {
        "peak"  : "1.1,  0,    2",
        "width" : "1.0,  0.01, 2",
        "skew"  : "-0.5, -2,   2",
      }
    else:
      raise ValueError(f"Unknown skewed Gaussian PDF subtype '{pdfSubType}'")
    defineSkewedGaussianPdf(
      fitManager, fitVariable, pdfName,
      parDefs,
      pdfSubType,
    )
  elif pdfTypeArgs[0] == "Histogram":
    defineHistogramPdf(
      fitManager, fitVariable, pdfName,
      {
        "smear" : "0" if "smear" in fixPars else "0,  0,    0.5",
        "shift" : "0" if "shift" in fixPars else "0, -0.25, 0.25",
        "scale" : "1" if "scale" in fixPars else "1,  0.5,  1.5",
      },
      outputDirName,
      templateDataFileName,
      templateDataTreeName,
      weightBranchName,
      comboIdName,
      cut =andCuts((cut, "(IsSignal == 0)")),
    )
  else:
    raise ValueError(f"Unknown background PDF type '{pdfTypeArgs[0]}'")

  bkgPdfWeightStartVal = 1.0
  fitManager.SetUp().LoadSpeciesPDF("BkgPdf", bkgPdfWeightStartVal)


def binnedTreeFilesIn(outputDirName: str) -> List[str]:
  """Returns list of file names with binned data"""
  binningFileName = f"{outputDirName}/DataBinsConfig.root"
  if not os.path.isfile(binningFileName):
    return []
  bins = ROOT.Bins("HSBins", binningFileName)
  binFileNames = [str(fileName) for fileName in bins.GetFileNames()]
  # verify that all bin files exist
  for binFileName in binFileNames:
    if not os.path.isfile(binFileName):
      return []
  return binFileNames


def setRooFitOptions(
  fitManager:       "ROOT.FitManager",
  nmbThreadsPerJob: int,
) -> None:
  """Sets general fit options"""
  print("Setting RooFit options")
  # see https://root.cern/doc/master/classRooAbsPdf.html#a52c4a5926a161bcb72eab46890b0590e
  fitManager.SetUp().AddFitOption(ROOT.RooFit.BatchMode(True))  # computes a batch of likelihood values at a time, uses faster math functions and possibly auto vectorization
                                                                # !Note! RooBatchCompute Library was revamped in ROOT 6.26/00 see https://github.com/root-project/root/tree/master/roofit/batchcompute
  fitManager.SetUp().AddFitOption(ROOT.RooFit.NumCPU(nmbThreadsPerJob))       # parallelizes calculation of likelihood using the given number of cores
  # fitManager.SetUp().AddFitOption(ROOT.RooFit.Parallelize(nmbThreadsPerJob))  # ROOT 6.28/00 global parallelization settings: enables use of RooFit's parallel minimization backend using the given number of cores
  fitManager.SetUp().AddFitOption(ROOT.RooFit.PrintLevel(2))
  # fitManager.SetUp().AddFitOption(ROOT.RooFit.Warnings(True))
  fitManager.SetUp().AddFitOption(ROOT.RooFit.Timer(True))  # times CPU and wall clock consumption of fit steps
  # force Minimizer
  # ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit")  # doesn't work
  # fitManager.SetUp().AddFitOption(ROOT.RooFit.Minimizer("Minuit"))  # overridden by HS::FIT::Minuit2::FitTo()
  # fitManager.SetMinimiser(ROOT.Minuit())
  # ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")
  # fitManager.SetUp().AddFitOption(ROOT.RooFit.Minimizer("Minuit2"))
  # fitManager.SetMinimiser(ROOT.Minuit2())
  # fitManager.SetUp().AddFitOption(ROOT.RooFit.Hesse(False))  # do not run HESSE
  # fitManager.SetUp().AddFitOption(ROOT.RooFit.Minos(True))  # run MINOS


def performFit(
  dataFileName:            str,                                            # path to file that contains data to fit
  outputDirName:           str,                                            # where to write all output files
  kinematicBinning:        Sequence[Tuple[str, int, float, float]],        # binning info with one tuple per dimension [ (<variable>, <nmb of bins>, <min value>, <max value>), ... ]
  pdfTypeSig:              str           = "Histogram",                    # type of signal PDF
  fixParsSig:              Sequence[str] = (),                             # fit-parameter names of signal function to fix
  pdfTypeBkg:              str           = "Histogram",                    # type of background PDF
  fixParsBkg:              Sequence[str] = (),                             # fit-parameter names of background function to fix
  commonCut:               str           = "",                             # optional selection cut(s) applied to data and template histograms
  dataCut:                 str           = "",                             # optional selection cut(s) applied to data only (in addition to commonCut)
  dataTreeName:            str           = "pippippimpimpmiss",            # name of tree that holds the data to fit
  templateDataSigFileName: str           = "",                             # name of file from which signal histogram is filled
  templateDataSigTreeName: str           = "pippippimpimpmiss",            # name of tree from which signal histogram is filled
  templateDataBkgFileName: str           = "",                             # name of file from which background histogram is filled
  templateDataBkgTreeName: str           = "pippippimpimpmiss",            # name of tree from which background histogram is filled
  fitVariable:             str           = "MissingMassSquared_Measured",  # name of branch that holds data to fit and template-data for signal and background, respectively
  fitRange:                str           = "-0.25, 3.75",                  # [(GeV/c)^2]
  comboIdName:             str           = "ComboID",                      # name of branch with unique combo ID
  regenBinnedTrees:        bool          = False,                          # if set, force regeneration of files with binned trees
  nmbThreadsPerJob:        int           = 5,
  nmbProofJobs:            int           = 20,  #TODO? automatically determine number of PROOF jobs
) -> None:
  """Sets up and performs fit"""
  # create the fit manager and set the output directory for fit results, plots, and weights
  fitDirName = None
  if not kinematicBinning:
    fitDirName = outputDirName
  else:
    # create subdirectory for binning
    binningVars = [binning[0] for binning in kinematicBinning]
    binningDirName = "_".join(binningVars)
    fitDirName = f"{outputDirName}/{binningDirName}"
  print(f"Fitting data in '{dataFileName}'" + ("" if not commonCut else f" applying cut '{commonCut}' to data and template histograms")
    + ("" if not dataCut else f" and additional cut '{dataCut}' to data")
    + (" using no binning" if not kinematicBinning else f" using binning '{kinematicBinning}'")
    + f" and writing output to '{fitDirName}'")
  gBenchmarkLabel = f"Time for fit in '{fitDirName}'"
  ROOT.gBenchmark.Start(gBenchmarkLabel)
  fitManager = ROOT.FitManager()
  fitManager.SetUp().SetOutDir(fitDirName)

  # define fit variable and set fit range
  print(f"Reading fit variable '{fitVariable}' from tree '{dataTreeName}' and using fit range {fitRange}")
  fitManager.SetUp().LoadVariable(f"{fitVariable}[{fitRange}]")
  # define combo-ID variable
  # the data tree must have a double branch of the given name containing a unique combo-ID number
  fitManager.SetUp().SetIDBranchName(comboIdName)

  # define kinematic bins
  # !Note! binning needs to be defined before any data are loaded
  for binning in kinematicBinning:
    fitManager.Bins().LoadBinVar(*binning)

  # define components of fit model
  if pdfTypeSig:
    defineSigPdf(
      fitManager, fitVariable, pdfTypeSig,
      fixPars              = fixParsSig,
      outputDirName        = fitDirName,
      templateDataFileName = templateDataSigFileName,
      templateDataTreeName = templateDataSigTreeName,
      weightBranchName     = "AccidWeightFactor",
      comboIdName          = comboIdName,
      cut                  = commonCut,
    )
  if pdfTypeBkg:
    defineBkgPdf(
      fitManager, fitVariable, pdfTypeBkg,
      fixPars              = fixParsBkg,
      outputDirName        = fitDirName,
      templateDataFileName = templateDataBkgFileName,
      templateDataTreeName = templateDataBkgTreeName,
      weightBranchName     = "AccidWeightFactor",
      comboIdName          = comboIdName,
      cut                  = commonCut
    )

  # create RF-sideband weights for data to be fitted
  rfSWeightLabel      = "RfSideband"
  rfSWeightObjectName = f"{rfSWeightLabel}WeightsData"
  rfSWeightFileName   = f"{fitDirName}/{rfSWeightObjectName}.root"
  readWeights(
    inputFileName     = dataFileName,
    inputTreeName     = dataTreeName,
    weightBranchName  = "AccidWeightFactor",
    comboIdName       = comboIdName,
    sWeightLabel      = rfSWeightLabel,
    sWeightFileName   = rfSWeightFileName,
    sWeightObjectName = rfSWeightObjectName,
    cut               = andCuts((commonCut, dataCut)),
  )
  # apply weights for RF-sideband subtraction
  fitManager.Data().LoadWeights(rfSWeightLabel, rfSWeightFileName, rfSWeightObjectName)

  # load and bin data to be fitted
  binFileNames = binnedTreeFilesIn(fitDirName) if kinematicBinning else []
  if not binFileNames or regenBinnedTrees:
    if kinematicBinning:
      if regenBinnedTrees:
        print("Forcing regeneration of binned tree files")
      else:
        print("Could not find (all) binned tree files; regenerating binned tree files")
    # fitManager.LoadData(dataTreeName, dataFileName, "Data")  # signature of HS::FIT::FitManager::LoadData() changed in b6e688a
    fitManager.LoadData(dataTreeName, dataFileName)
  else:
    print("Using existing binned tree files:")
    for binFileName in binFileNames:
      print(f"    {binFileName}")
    fitManager.ReloadData(dataTreeName, dataFileName, "Data")

  # perform fit and create fit-result plots
  if not kinematicBinning:
    nmbThreadsPerJob = 5 * nmbThreadsPerJob
  setRooFitOptions(fitManager, nmbThreadsPerJob)
  print("Using the following global fit options")
  fitManager.SetUp().FitOptions().Print("")
  if kinematicBinning:
    fitManager.SetRedirectOutput()  # redirect console output to files
    print(f"Running {nmbProofJobs} PROOF jobs each with {nmbThreadsPerJob} threads in parallel")
    ROOT.gEnv.SetValue("ProofLite.Sandbox", "$PWD/.proof/")
    ROOT.Proof.Go(fitManager, nmbProofJobs)
  else:
    print(f"Running {nmbThreadsPerJob} threads in parallel")
    ROOT.Here.Go(fitManager)

  fitManager.WriteThis()  # write to disk
  ROOT.gBenchmark.Show(gBenchmarkLabel)


if __name__ == "__main__":
  plotTools.printGitInfo()
  os.nice(18)  # run all processes with second highest niceness level
  ROOT.gROOT.SetBatch(True)
  plotTools.setupPlotStyle()
  ROOT.gROOT.ProcessLine(f".x {os.environ['BRUFIT']}/macros/LoadBru.C")
  ROOT.gBenchmark.Start("Total execution time")

  # echo and parse command line
  # bggenFileName = f"./data/MCbggen/2017_01-ver03/pippippimpimpmiss_flatTree.MCbggen_2017_01-ver03.root.brufit"
  bggenFileName = f"./data/MCbggen/2018_01-ver02/pippippimpimpmiss_flatTree.MCbggen_2018_01-ver02.root.brufit"
  # dataFileName  = bggenFileName
  # dataFileName  = f"./data/RD/2017_01-ver04/pippippimpimpmiss_flatTree.RD_2017_01-ver04_030730.root.brufit"
  dataFileName  = f"./data/RD/2018_01-ver02/old/pippippimpimpmiss_flatTree.RD_2018_01-ver02_041003.root.brufit"
  # dataFileName  = f"./data/RD/2018_01-ver02/old/pippippimpimpmiss_flatTree.RD_2018_01-ver02_042030.root.brufit"
  # dataFileName  = f"./data/RD/2018_01-ver02/old/pippippimpimpmiss_flatTree.RD_2018_01-ver02_042550.root.brufit"
  print(f"Script was called using: '{' '.join(sys.argv)}'")
  parser = argparse.ArgumentParser(description="Plots BruFit results.")
  parser.add_argument("outputDirName", nargs = "?", type = str, default = "./BruFitOutput", help = "The path to the BruFit output directory; (default: '%(default)s')")
  parser.add_argument("dataFileName",  nargs = "?", type = str, default = dataFileName,     help = "The path to the input data tree in BruFit format; (default: '%(default)s')")
  parser.add_argument("bggenFileName", nargs = "?", type = str, default = bggenFileName,    help = "The path to the input bggen MC tree in BruFit format; (default: '%(default)s')")
  parser.add_argument("--pdfTypeSig",               type = str, default = "Histogram",      help = "Type of signal PDF to use in fit; (default: '%(default)s')")
  parser.add_argument("--fixParsSig",  nargs = "*", type = str, default = [],               help = "Names of parameters of signal PDF to fix in fit; (default: none)")
  parser.add_argument("--pdfTypeBkg",               type = str, default = "Histogram",      help = "Type of background PDF to use in fit; (default: '%(default)s')")
  parser.add_argument("--fixParsBkg",  nargs = "*", type = str, default = [],               help = "Names of parameters of background PDF to fix in fit; (default: none)")
  args = parser.parse_args()
  dataCut           = ""
  # dataCut           = "(IsSignal == 1)"  # fit bggen signal data
  # dataCut           = "(IsSignal == 0)"  # fit bggen background data
  # additionalCut     = ""
  # additionalCut     = "(NmbUnusedShowers == 0)"
  additionalCut     = "(NmbUnusedShowers == 0) && (MissingProtonP > 0.5)"
  kinematicBinnings: List[List[Tuple[str, int, float, float]]] = [
    [],  # no binning -> fit overall distribution
    # # 1D binnings
    # [("BeamEnergy",          90,    2.9,   11.9)],  # [GeV]
    # [("MissingProtonP",     100,    0,      5)],    # [GeV/c]
    # [("MissingProtonTheta",  72,    0,     90)],    # [deg]
    # [("MissingProtonPhi",    72, -180,   +180)],    # [deg]
    # 2D binnings
    # [
    #   ("MissingProtonTheta", 9, 0, 90),  # [deg]
    #   ("MissingProtonP",    25, 0,  5),  # [GeV/c]
    # ],
    [
      ("MissingProtonTheta", 10, 0, 20),  # [deg]
      ("MissingProtonP",      9, 0,  9),  # [GeV/c]
    ],
  ]

  dataSets = {
    "Total"   : "",
    "Found"   : "(TrackFound == 1)",
    "Missing" : "(TrackFound == 0)",
  }

  #TODO calculate chi^2
  # see https://root.cern/doc/master/rf109__chi2residpull_8py.html
  # and https://root-forum.cern.ch/t/how-to-correctly-extract-the-chi2-ndf-p-value-of-a-roofitresult/45956
  # fit all datasets and bins
  for dataSetName, dataSetCut in dataSets.items():
    if kinematicBinnings:
      for kinematicBinning in kinematicBinnings:
        performFit(
          args.dataFileName,
          f"{args.outputDirName}/{dataSetName}",
          kinematicBinning,
          pdfTypeSig              = args.pdfTypeSig,
          fixParsSig              = args.fixParsSig,
          pdfTypeBkg              = args.pdfTypeBkg,
          fixParsBkg              = args.fixParsBkg,
          commonCut               = andCuts((dataSetCut, additionalCut)),
          dataCut                 = dataCut,
          templateDataSigFileName = args.bggenFileName,
          templateDataBkgFileName = args.bggenFileName,
        )

  ROOT.gBenchmark.Show("Total execution time")
