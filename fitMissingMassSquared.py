#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import os

import ROOT

import makePlots  # defines helper functions to generate histograms from data trees


#TODO convert functions to class with member functions
def readWeights(
  inputFileName,
  inputTreeName,
  weightBranchName,
  comboIdName,
  sWeightLabel,       # label used when applying weights
  sWeightFileName,    # ROOT file name to which weight object will be written
  sWeightObjectName,  # name of weight object saved to ROOT file
  cut = None          # optional selection cut to apply
):
  print(f"Reading weights '{weightBranchName}' from tree '{inputTreeName}' in file '{inputFileName}'"
  f", writing them to key '{sWeightObjectName}' in file '{sWeightFileName}', and assigning label '{sWeightLabel}'"
  + ("" if cut is None else f" while applying cut(s) '{cut}'"), flush = True)
  currentDir = ROOT.gDirectory
  inputFile = ROOT.TFile.Open(inputFileName, "READ")
  inputTree = inputFile.Get(inputTreeName)
  currentDir.cd()
  weights = ROOT.Weights(sWeightObjectName)  # name of the Weights object
  weights.SetFile(sWeightFileName)
  weights.SetSpecies(sWeightLabel)
  weights.SetIDName(comboIdName)
  weights.WeightBySelection(inputTree, ("(1)" if cut is None else cut), weightBranchName)
  weights.SortWeights()
  weights.Save()


def defineSigPdf(
  fitManager,
  fitVariable,
  outputDirName        = None,
  templateDataFileName = None,
  templateDataTreeName = None,
  weightBranchName     = None,
  comboIdName          = None,
  cut                  = None
):
  #TODO implement function that constructs a selectable PDF and can be used construct signal as well as background PDFs
  print("Defining signal PDF 'SigPdf'", flush = True)

  # # single Gaussian
  # fitManager.SetUp().FactoryPDF(f"Gaussian::SigPdf({fitVariable}, mean_SigPdf[{meanStartVal}, 0, 2], width_SigPdf[{widthStartVal}, 0.1, 3])")

  # # double Gaussian
  # # with same mean
  # fitManager.SetUp().FactoryPDF("SUM::SigPdf("
  #   f"r_SigPdf[0.5, 0, 1] * Gaussian::SigPdf_N1({fitVariable}, mean_SigPdf[{0.9383**2}, 0, 2], width1_SigPdf[1.0, 0.01, 2]),"  # wide Gaussian
  #                         f"Gaussian::SigPdf_N2({fitVariable}, mean_SigPdf,                    width2_SigPdf[0.2, 0.01, 2])"   # narrow Gaussian
  # ")")
  # # with separate means
  # fitManager.SetUp().FactoryPDF("SUM::SigPdf("
  #   f"r_SigPdf[0.5, 0, 1] * Gaussian::SigPdf_N1({fitVariable}, mean1_SigPdf[1.0,         0, 2], width1_SigPdf[1.0, 0.01, 2]),"  # wide Gaussian
  #                         f"Gaussian::SigPdf_N2({fitVariable}, mean2_SigPdf[{0.9383**2}, 0, 2], width2_SigPdf[0.2, 0.01, 2])"   # narrow Gaussian
  # ")")

  # create RF-sideband weights for histogram PDF
  rfSWeightLabel      = "RfSideband"
  rfSWeightFileName   = f"{outputDirName}/rfSidebandWeightsSigPdf.root"
  rfSWeightObjectName = f"{rfSWeightLabel}WeightsSigPdf"
  readWeights(
    inputFileName     = templateDataFileName,
    inputTreeName     = templateDataTreeName,
    weightBranchName  = weightBranchName,
    comboIdName       = comboIdName,
    sWeightLabel      = rfSWeightLabel,
    sWeightFileName   = rfSWeightFileName,
    sWeightObjectName = rfSWeightObjectName,
    # cut               = cut
    cut               = ("" if cut is None else f"({cut}) && ") + '(IsSignal == 1)'
  )
  # histogram PDF with fudge parameters that allow small deviations from the given shape:
  #     smear  = width of Gaussian the histogram is convoluted with
  #     offset = shift of histogram in x-direction
  #     scale  = scaling factor of histogram in x-direction
  fitManager.SetUp().FactoryPDF(
    f"RooHSEventsHistPDF::SigPdf({fitVariable}, smear_SigPdf[0, 0, 0.5], offset_SigPdf[0, -0.25, 0.25], scale_SigPdf[1, 0.5, 1.5])"  # correct
    # f"RooHSEventsHistPDF::SigPdf({fitVariable}, smear_SigPdf[0], offset_SigPdf[0, -0.25, 0.25], scale_SigPdf[1, 0.5, 1.5])"  # line at average
    # f"RooHSEventsHistPDF::SigPdf({fitVariable}, smear_SigPdf[0, 0, 0.5], offset_SigPdf[0], scale_SigPdf[1, 0.5, 1.5])"  # line at 0
    # f"RooHSEventsHistPDF::SigPdf({fitVariable}, smear_SigPdf[0, 0, 0.5], offset_SigPdf[0, -0.25, 0.25], scale_SigPdf[1])"  # weird
    # f"RooHSEventsHistPDF::SigPdf({fitVariable}, smear_SigPdf[0, 0, 0.5], offset_SigPdf[0], scale_SigPdf[1])"  # line
    # f"RooHSEventsHistPDF::SigPdf({fitVariable}, smear_SigPdf[0], offset_SigPdf[0, -0.25, 0.25], scale_SigPdf[1])"  # line above peak
    # f"RooHSEventsHistPDF::SigPdf({fitVariable}, smear_SigPdf[0], offset_SigPdf[0], scale_SigPdf[1, 0.5, 1.5])"  # line
    # f"RooHSEventsHistPDF::SigPdf({fitVariable}, smear_SigPdf[0], offset_SigPdf[0], scale_SigPdf[1])"  # line at 0
    # f"RooHSEventsHistPDF::SigPdf({fitVariable}, smear_SigPdf[0], offset_SigPdf[0.25], scale_SigPdf[1])"  # line at 0
    + f"WEIGHTS@{rfSWeightLabel},{rfSWeightFileName},{rfSWeightObjectName}"  # apply sWeights created above; !Note! no whitespace allowed in this string
  )
  # constrain PDF fudge parameters
  # the constraints are derived from the smear, offset, and scale parameter initial values and limits, i.e.
  #     Gaussian mean  = initial value
  #     Gaussian sigma = max / 5 for smear
  #     Gaussian sigma = range / 10 for offset and scale
  sigPdf = fitManager.SetUp().WS().pdf("SigPdf")
  fitManager.SetUp().AddGausConstraint(sigPdf.AlphaConstraint())  # constrain smear
  fitManager.SetUp().AddGausConstraint(sigPdf.OffConstraint())    # constrain offset
  fitManager.SetUp().AddGausConstraint(sigPdf.ScaleConstraint())  # constrain scale
  # load data for template histograms
  print(f"Loading data for histogram PDF 'SigPdf' from tree '{templateDataTreeName}' in file '{templateDataFileName}'"
  + f" and applying weights from '{rfSWeightObjectName}' in file '{rfSWeightFileName}'", flush = True)
  fitManager.LoadSimulated(templateDataTreeName, templateDataFileName, "SigPdf")

  sigPdfWeightStartVal = 1.0
  fitManager.SetUp().LoadSpeciesPDF("SigPdf", sigPdfWeightStartVal)


def defineBkgPdf(
  fitManager,
  fitVariable,
  outputDirName        = None,
  templateDataFileName = None,
  templateDataTreeName = None,
  weightBranchName     = None,
  comboIdName          = None,
  cut                  = None
):
  print("Defining background PDF 'BkgPdf'", flush = True)

  # # 2nd-order positive-definite polynomial
  # fitManager.SetUp().FactoryPDF(f"GenericPdf::BkgPdf('@1 * @1 + (@2 + @3 * @0) * (@2 + @3 * @0)', {{{fitVariable}, p0_BkgPdf[0, -100, 100], p1_BkgPdf[0, -100, 100], p2_BkgPdf[0, -100, 100]}})")

  # Chebychev polynomial
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BkgPdf({fitVariable}, {{}})")
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BkgPdf({fitVariable}, {{p0_BkgPdf[0, 0, 100]}})")
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BkgPdf({fitVariable}, {{p0_BkgPdf[0, -1, 1], p1_BkgPdf[0, -1, 1]}})")
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BkgPdf({fitVariable}, {{p0_BkgPdf[0, -1, 1], p1_BkgPdf[0, -1, 1], p2_BkgPdf[0, -1, 1]}})")

  # # Bernstein polynomial
  # # see https://root.cern.ch/doc/master/classRooBernstein.html
  # fitManager.SetUp().FactoryPDF(f"Bernstein::BkgPdf({fitVariable}, {{}})")
  # fitManager.SetUp().FactoryPDF(f"Bernstein::BkgPdf({fitVariable}, {{p0_BkgPdf[0, 0, 1]}})")
  # fitManager.SetUp().FactoryPDF(f"Bernstein::BkgPdf({fitVariable}, {{p0_BkgPdf[0, 0, 1], p1_BkgPdf[0, 0, 1]}})")
  # fitManager.SetUp().FactoryPDF(f"Bernstein::BkgPdf({fitVariable}, {{p0_BkgPdf[0, 0, 1], p1_BkgPdf[0, 0, 1], p2_BkgPdf[0, 0, 1]}})")

  # # double Gaussian with separate means
  # fitManager.SetUp().FactoryPDF("SUM::BkgPdf("
  #   f"r_BkgPdf[0.5, 0, 1] * Gaussian::BkgPdf_N1({fitVariable}, mean1_BkgPdf[1.1, 0, 2], width1_BkgPdf[1.0, 0.01, 2]),"  # wide Gaussian
  #                         f"Gaussian::BkgPdf_N2({fitVariable}, mean2_BkgPdf[0.9, 0, 2], width2_BkgPdf[0.5, 0.01, 2])"   # narrow Gaussian
  # ")")

  # various implementations of skewed Gaussians
  # # "logarithmic Gaussian"
  # # Eq. (9) in NIMA 441 (2000) 401 at https://doi.org/10.1016/S0168-9002(99)00992-4
  # #TODO does not work; RooFit complains about PDF being -inf; unclear why
  # fitManager.SetUp().FactoryPDF(f"Novosibirsk::BkgPdf({fitVariable}, peak_BkgPdf[1.1, 0, 2], width_BkgPdf[1.0, 0.01, 2], skew_BkgPdf[-0.5, -2, 2])")
  # # skew normal PDF (see https://en.wikipedia.org/wiki/Skew_normal_distribution)
  # fitManager.SetUp().FactoryPDF("EXPR::BkgPdf("
  #   f"'2 * TMath::Gaus({fitVariable}, peak_BkgPdf, width_BkgPdf, true)"
  #   f"* TMath::Erfc(skew_BkgPdf * (({fitVariable} - peak_BkgPdf) / width_BkgPdf))',"
  #   f"{fitVariable}, peak_BkgPdf[0.8, 0, 2], width_BkgPdf[0.9, 0.01, 2], skew_BkgPdf[-2.5, -5, 5])")
  # # exponentially modified Gaussian PDF (see https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution)
  # # smaller skew means higher skewness
  # # max = peak + 1 / skew
  # fitManager.SetUp().FactoryPDF("EXPR::BkgPdf("
  #   f"'(skew_BkgPdf / 2) * TMath::Exp((skew_BkgPdf / 2) * (2 * peak_BkgPdf + skew_BkgPdf * width_BkgPdf * width_BkgPdf - 2 * {fitVariable}))"
  #   f"* TMath::Erfc((peak_BkgPdf + skew_BkgPdf * width_BkgPdf * width_BkgPdf - {fitVariable}) / (sqrt(2) * width_BkgPdf))',"
  #   f"{fitVariable}, peak_BkgPdf[0.8, 0, 2], width_BkgPdf[0.5, 0.01, 2], skew_BkgPdf[2, 0, 15])")
  #TODO try https://en.wikipedia.org/wiki/Normal-exponential-gamma_distribution

  # #TODO move code to define histogram PDF to separate function
  # create RF-sideband weights for histogram PDF
  rfSWeightLabel      = "RfSideband"
  rfSWeightFileName   = f"{outputDirName}/rfSidebandWeightsBkgPdf.root"
  rfSWeightObjectName = f"{rfSWeightLabel}WeightsBkgPdf"
  readWeights(
    inputFileName     = templateDataFileName,
    inputTreeName     = templateDataTreeName,
    weightBranchName  = weightBranchName,
    comboIdName       = comboIdName,
    sWeightLabel      = rfSWeightLabel,
    sWeightFileName   = rfSWeightFileName,
    sWeightObjectName = rfSWeightObjectName,
    # cut               = cut
    cut               = ("" if cut is None else f"({cut}) && ") + '(IsSignal == 0)'
  )
  # histogram PDF with fudge parameters that allow small deviations from the given shape:
  #     smear  = width of Gaussian the histogram is convoluted with
  #     offset = shift of histogram in x-direction
  #     scale  = scaling factor of histogram in x-direction
  fitManager.SetUp().FactoryPDF(
    f"RooHSEventsHistPDF::BkgPdf({fitVariable}, smear_BkgPdf[0, 0, 0.5], offset_BkgPdf[0, -0.25, 0.25], scale_BkgPdf[1, 0.5, 1.5])"
    # f"RooHSEventsHistPDF::BkgPdf({fitVariable}, smear_BkgPdf[0], offset_BkgPdf[0, -0.25, 0.25], scale_BkgPdf[1, 0.5, 1.5])"
    + f"WEIGHTS@{rfSWeightLabel},{rfSWeightFileName},{rfSWeightObjectName}"  # apply sWeights created above; !Note! no whitespace allowed in this string
  )
  # constrain PDF fudge parameters
  # the constraints are derived from the smear, offset, and scale parameter initial values and limits, i.e.
  #     Gaussian mean  = initial value
  #     Gaussian sigma = max / 5 for smear
  #     Gaussian sigma = range / 10 for offset and scale
  bkgPdf = fitManager.SetUp().WS().pdf("BkgPdf")
  fitManager.SetUp().AddGausConstraint(bkgPdf.AlphaConstraint())  # constrain smear
  fitManager.SetUp().AddGausConstraint(bkgPdf.OffConstraint())    # constrain offset
  fitManager.SetUp().AddGausConstraint(bkgPdf.ScaleConstraint())  # constrain scale
  # load data for template histograms
  print(f"Loading data for histogram PDF 'BkgPdf' from tree '{templateDataTreeName}' in file '{templateDataFileName}'"
  + f" and applying weights from '{rfSWeightObjectName}' in file '{rfSWeightFileName}'", flush = True)
  fitManager.LoadSimulated(templateDataTreeName, templateDataFileName, "BkgPdf")

  bkgPdfWeightStartVal = 1.0
  fitManager.SetUp().LoadSpeciesPDF("BkgPdf", bkgPdfWeightStartVal)


def binnedTreeFilesIn(outputDirName):
  binningFileName = f"{outputDirName}/DataBinsConfig.root"
  if not os.path.isfile(binningFileName):
    return None
  bins = ROOT.Bins("HSBins", binningFileName)
  binFileNames = [str(fileName) for fileName in bins.GetFileNames()]
  # verify that all bin files exist
  for binFileName in binFileNames:
    if not os.path.isfile(binFileName):
      return None
  return binFileNames


def setRooFitOptions(
  fitManager,
  nmbThreadsPerJob
):
  print("Setting RooFit options", flush = True)
  # see https://root.cern/doc/master/classRooAbsPdf.html#a52c4a5926a161bcb72eab46890b0590e
  fitManager.SetUp().AddFitOption(ROOT.RooFit.BatchMode(True))  # computes a batch of likelihood values at a time, uses faster math functions and possibly auto vectorization
                                                                # !Note! RooBatchCompute Library was revamped in ROOT 6.26/00 see https://github.com/root-project/root/tree/master/roofit/batchcompute
  fitManager.SetUp().AddFitOption(ROOT.RooFit.NumCPU(nmbThreadsPerJob))       # parallelizes calculation of likelihood using the given number of cores
  # fitManager.SetUp().AddFitOption(ROOT.RooFit.Parallelize(nmbThreadsPerJob))  # ROOT 6.28/00 global parallelization settings: enables use of RooFit's parallel minimization backend using the given number of cores
  fitManager.SetUp().AddFitOption(ROOT.RooFit.PrintLevel(2))
  fitManager.SetUp().AddFitOption(ROOT.RooFit.Timer(True))  # times CPU and wall clock consumption of fit steps
  # fitManager.SetUp().AddFitOption(ROOT.RooFit.Hesse(False))  # do not run HESSE
  # fitManager.SetUp().AddFitOption(ROOT.RooFit.Minos(True))  # run MINOS


def performFit(
  dataFileName,      # path to data to fit
  outputDirName,     # where to write all output files
  kinematicBinning,  # list of with binning info, one tuple per dimension [ (<variable>, <nmb of bins>, <min value>, <max value>), ... ]
  commonCut               = None,                 # optional selection cut(s) applied to data and template histograms
  dataCut                 = None,                 # optional selection cut(s) applied to data only (in addition to commonCut)
  dataTreeName            = "pippippimpimpmiss",  # tree name of data to fit
  templateDataSigFileName = None,                 # path to template data for signal histogram
  templateDataSigTreeName = "pippippimpimpmiss",  # tree name of data for signal histogram template
  templateDataBkgFileName = None,                 # path to template data for background histogram
  templateDataBkgTreeName = "pippippimpimpmiss",  # tree name of data for background histogram template
  fitVariable             = "MissingMassSquared_Measured",  # branch name that contains data to fit and template-data for signal and background, respectively
  fitRange                = "-0.25, 3.75",  # [(GeV/c)^2]
  comboIdName             = "ComboID",  # branch name with unique ID for each combo
  regenBinnedTrees        = False,  # if set, force regeneration of files with binned trees
  nmbThreadsPerJob        = 5,
  nmbProofJobs            = 10  #TODO? automatically determine number of PROOF jobs
):
  # create the fit manager and set the output directory for fit results, plots, and weights
  fitDirName = None
  if kinematicBinning is None:
    fitDirName = outputDirName
  else:
    # create subdirectory for binning
    binningVars = [binning[0] for binning in kinematicBinning]
    binningDirName = "_".join(binningVars)
    fitDirName = f"{outputDirName}/{binningDirName}"
  print(f"fitting data in '{dataFileName}'" + ("" if commonCut is None else f" applying cut '{commonCut}' to data and template histograms")
    + ("" if dataCut is None else f" and additional cut '{dataCut}' to data")
    + (" using no binning" if kinematicBinning is None else f" using binning '{kinematicBinning}'")
    + f" and writing output to '{fitDirName}'")
  gBenchmarkLabel = f"Time for fit in '{fitDirName}'"
  ROOT.gBenchmark.Start(gBenchmarkLabel)
  fitManager = ROOT.FitManager()
  fitManager.SetUp().SetOutDir(fitDirName)

  # define fit variable and set fit range
  print(f"reading fit variable '{fitVariable}' from tree '{dataTreeName}' and using fit range {fitRange}")
  fitManager.SetUp().LoadVariable(f"{fitVariable}[{fitRange}]")
  # define combo-ID variable
  # the data tree must have a double branch of the given name containing a unique combo-ID number
  fitManager.SetUp().SetIDBranchName(comboIdName)

  # define kinematic bins
  # !Note! binning needs to be defined before any data are loaded
  if kinematicBinning:
    for binning in kinematicBinning:
      fitManager.Bins().LoadBinVar(*binning)

  # define components of fit model
  defineSigPdf(
    fitManager, fitVariable,
    outputDirName        = fitDirName,
    templateDataFileName = templateDataSigFileName,
    templateDataTreeName = templateDataSigTreeName,
    weightBranchName     = "AccidWeightFactor",
    comboIdName          = comboIdName,
    cut                  = commonCut
  )
  defineBkgPdf(
    fitManager, fitVariable,
    outputDirName        = fitDirName,
    templateDataFileName = templateDataBkgFileName,
    templateDataTreeName = templateDataBkgTreeName,
    weightBranchName     = "AccidWeightFactor",
    comboIdName          = comboIdName,
    cut                  = commonCut
  )

  # create RF-sideband weights for data to be fitted
  rfSWeightLabel      = "RfSideband"
  rfSWeightFileName   = f"{fitDirName}/rfSidebandWeightsData.root"
  rfSWeightObjectName = f"{rfSWeightLabel}WeightsData"
  readWeights(
    inputFileName     = dataFileName,
    inputTreeName     = dataTreeName,
    weightBranchName  = "AccidWeightFactor",
    comboIdName       = comboIdName,
    sWeightLabel      = rfSWeightLabel,
    sWeightFileName   = rfSWeightFileName,
    sWeightObjectName = rfSWeightObjectName,
    cut               = ("" if commonCut is None else f"({commonCut}) && ") + f"({dataCut})"
  )
  # apply weights for RF-sideband subtraction
  fitManager.Data().LoadWeights(rfSWeightLabel, rfSWeightFileName, rfSWeightObjectName)

  # load and bin data to be fitted
  binFileNames = binnedTreeFilesIn(fitDirName) if kinematicBinning else None
  if binFileNames is None or regenBinnedTrees:
    if kinematicBinning:
      if regenBinnedTrees:
        print("Forcing regeneration of binned tree files", flush = True)
      else:
        print("Could not find (all) binned tree files; regenerating binned tree files", flush = True)
    fitManager.LoadData(dataTreeName, dataFileName, "Data")
  else:
    print("Using existing binned tree files:", flush = True)
    for binFileName in binFileNames:
      print(f"    {binFileName}", flush = True)
    fitManager.ReloadData(dataTreeName, dataFileName, "Data")

  # perform fit and create fit-result plots
  if not kinematicBinning:
    nmbThreadsPerJob = 8 * nmbThreadsPerJob
  setRooFitOptions(fitManager, nmbThreadsPerJob)
  if kinematicBinning:
    fitManager.SetRedirectOutput()  # redirect console output to files
    print(f"running {nmbProofJobs} PROOF jobs each with {nmbThreadsPerJob} threads in parallel", flush = True)
    ROOT.Proof.Go(fitManager, nmbProofJobs)
  else:
    print(f"running {nmbThreadsPerJob} threads in parallel", flush = True)
    ROOT.Here.Go(fitManager)

  fitManager.WriteThis()  # write to disk
  ROOT.gBenchmark.Show(gBenchmarkLabel)


if __name__ == "__main__":
  os.nice(18)  # run all processes with second highest niceness level
  ROOT.gROOT.SetBatch(True)
  makePlots.setupPlotStyle()
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable
  ROOT.gBenchmark.Start("Total execution time")

  # dataSample        = "030730",
  dataSample        = "bggen_2017_01-ver03"
  dataFileName      = f"../ReactionEfficiency/pippippimpimpmiss_flatTree.{dataSample}.root.brufit"
  # dataCut           = None
  dataCut           = "(IsSignal == 1)"  # fit bggen signal data
  # dataCut           = "(IsSignal == 0)"  # fit bggen background data
  outputDirName     = "./BruFitOutput"
  kinematicBinnings = [
    None,  # no binning -> fit overall distribution
    # 1D binnings; only one binning par variable name allowed
    [("BeamEnergy",          9,    3.0,   12.0)],
    [("MissingProtonP",     10,    0.0,    3.5)],
    [("MissingProtonTheta", 10,    0.0,   65.0)],
    [("MissingProtonPhi",   10, -180.0, +180.0)]
  ]

  dataSets = {
    "Total"   : None,
    "Found"   : "(TrackFound == 1)",
    "Missing" : "(TrackFound == 0)"
  }

  #TODO calculate chi^2
  # see https://root.cern/doc/master/rf109__chi2residpull_8py.html
  # and https://root-forum.cern.ch/t/how-to-correctly-extract-the-chi2-ndf-p-value-of-a-roofitresult/45956
  # fit all datasets and bins
  for dataSetName, dataSetCut in dataSets.items():
    if kinematicBinnings:
      for kinematicBinning in kinematicBinnings:
        performFit(
          dataFileName,
          f"{outputDirName}/{dataSetName}",
          kinematicBinning,
          commonCut               = dataSetCut,
          dataCut                 = dataCut,
          templateDataSigFileName = dataFileName,
          templateDataBkgFileName = dataFileName
        )

  ROOT.gBenchmark.Show("Total execution time")
