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
  #TODO make different PDFs selectable
  print("Defining signal PDF 'SigPdf'", flush = True)

  # # define single Gaussian
  # fitManager.SetUp().FactoryPDF(f"Gaussian::SigPdf({fitVariable}, mean_SigPdf[{meanStartVal}, 0, 2], width_SigPdf[{widthStartVal}, 0.1, 3])")

  # # define double Gaussian
  # # separate means
  # fitManager.SetUp().FactoryPDF("SUM::SigPdf("
  #   f"r_SigPdf[0.5, 0, 1] * Gaussian::SigPdf_N1({fitVariable}, mean1_SigPdf[1.0,         0, 2], width1_SigPdf[1.0, 0.01, 2]),"  # wide Gaussian
  #                         f"Gaussian::SigPdf_N2({fitVariable}, mean2_SigPdf[{0.9383**2}, 0, 2], width2_SigPdf[0.2, 0.01, 2])"   # narrow Gaussian
  #   ")")
  # # same mean
  # fitManager.SetUp().FactoryPDF("SUM::SigPdf("
  #   f"r_SigPdf[0.5, 0, 1] * Gaussian::SigPdf_N1({fitVariable}, mean_SigPdf[{0.9383**2}, 0, 2], width1_SigPdf[1.0, 0.01, 2]),"  # wide Gaussian
  #                         f"Gaussian::SigPdf_N2({fitVariable}, mean_SigPdf,                    width2_SigPdf[0.2, 0.01, 2])"   # narrow Gaussian
  #   ")")

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
  # define histogram PDF with fudge parameters that allow small deviations from the given shape:
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
  #TODO make different PDFs selectable
  print("Defining background PDF 'BkgPdf'", flush = True)

  # # define 2nd-order positive-definite polynomial
  # fitManager.SetUp().FactoryPDF(f"GenericPdf::BkgPdf('@1 * @1 + (@2 + @3 * @0) * (@2 + @3 * @0)', {{{fitVariable}, p0_BkgPdf[0, -100, 100], p1_BkgPdf[0, -100, 100], p2_BkgPdf[0, -100, 100]}})")

  # define Chebychev polynomial
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BkgPdf({fitVariable}, {{}})")
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BkgPdf({fitVariable}, {{p0_BkgPdf[0, 0, 100]}})")
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BkgPdf({fitVariable}, {{p0_BkgPdf[0, -1, 1], p1_BkgPdf[0, -1, 1]}})")
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BkgPdf({fitVariable}, {{p0_BkgPdf[0, -1, 1], p1_BkgPdf[0, -1, 1], p2_BkgPdf[0, -1, 1]}})")

  # # define Bernstein polynomial
  # # see https://root.cern.ch/doc/master/classRooBernstein.html
  # fitManager.SetUp().FactoryPDF(f"Bernstein::BkgPdf({fitVariable}, {{}})")
  # fitManager.SetUp().FactoryPDF(f"Bernstein::BkgPdf({fitVariable}, {{p0_BkgPdf[0, 0, 1]}})")
  # fitManager.SetUp().FactoryPDF(f"Bernstein::BkgPdf({fitVariable}, {{p0_BkgPdf[0, 0, 1], p1_BkgPdf[0, 0, 1]}})")
  # fitManager.SetUp().FactoryPDF(f"Bernstein::BkgPdf({fitVariable}, {{p0_BkgPdf[0, 0, 1], p1_BkgPdf[0, 0, 1], p2_BkgPdf[0, 0, 1]}})")

  #TODO move code to define histogram PDF to separate function
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
  # define histogram PDF with fudge parameters that allow small deviations from the given shape:
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
  dataFileName,       # path to data to fit
  outputDirName,      # where to write all output files
  kinematicBinnings,  # list of tuples with binning info [ (<variable>, <nmb of bins>, <min value>, <max value>), ... ]
  cut                     = None,                 # optional selection cut(s) to apply to data and template histograms
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
  nmbProofJobs            = 9  #TODO? automatically determine number of PROOF jobs
):
  print(f"fitting data in '{dataFileName}'" + ("" if cut is None else f" applying cut '{cut}'")
    + (" using no binning" if kinematicBinnings is None else f" using binning '{kinematicBinnings}'")
    + f" and writing output to '{outputDirName}'")
  print(f"reading fit variable '{fitVariable}' from tree '{dataTreeName}' and using fit range {fitRange}")
  gBenchmarkLabel = f"Time for fit in '{outputDirName}'"
  ROOT.gBenchmark.Start(gBenchmarkLabel)

  # create the fit manager and set the output directory for fit results, plots, and weights
  fitManager = ROOT.FitManager()
  fitManager.SetUp().SetOutDir(outputDirName)
  # define fit variable and set fit range
  fitManager.SetUp().LoadVariable(f"{fitVariable}[{fitRange}]")
  # define combo-ID variable
  # the data tree must have a double branch of the given name containing a unique combo-ID number
  fitManager.SetUp().SetIDBranchName(comboIdName)

  # define kinematic bins
  # !Note! binning needs to be defined before any data are loaded
  if kinematicBinnings:
    for binning in kinematicBinnings:
      fitManager.Bins().LoadBinVar(*binning)

  # define components of fit model
  defineSigPdf(
    fitManager, fitVariable,
    outputDirName        = outputDirName,
    templateDataFileName = templateDataSigFileName,
    templateDataTreeName = templateDataSigTreeName,
    weightBranchName     = "AccidWeightFactor",
    comboIdName          = comboIdName,
    cut                  = cut
  )
  # defineBkgPdf(
  #   fitManager, fitVariable,
  #   outputDirName        = outputDirName,
  #   templateDataFileName = templateDataBkgFileName,
  #   templateDataTreeName = templateDataBkgTreeName,
  #   weightBranchName     = "AccidWeightFactor",
  #   comboIdName          = comboIdName,
  #   cut                  = cut
  # )

  # create RF-sideband weights for data to be fitted
  rfSWeightLabel      = "RfSideband"
  rfSWeightFileName   = f"{outputDirName}/rfSidebandWeightsData.root"
  rfSWeightObjectName = f"{rfSWeightLabel}WeightsData"
  readWeights(
    inputFileName     = dataFileName,
    inputTreeName     = dataTreeName,
    weightBranchName  = "AccidWeightFactor",
    comboIdName       = comboIdName,
    sWeightLabel      = rfSWeightLabel,
    sWeightFileName   = rfSWeightFileName,
    sWeightObjectName = rfSWeightObjectName,
    cut               = cut
    # cut               = ("" if cut is None else f"({cut}) && ") + '(IsSignal == 0)'
  )
  # apply weights for RF-sideband subtraction
  fitManager.Data().LoadWeights(rfSWeightLabel, rfSWeightFileName, rfSWeightObjectName)

  # load and bin data to be fitted
  binFileNames = binnedTreeFilesIn(outputDirName) if kinematicBinnings else None
  if binFileNames is None or regenBinnedTrees:
    if kinematicBinnings:
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

  # perform fit and plot fit result
  if not kinematicBinnings:
    nmbThreadsPerJob = 8 * nmbThreadsPerJob
  setRooFitOptions(fitManager, nmbThreadsPerJob)
  if kinematicBinnings:
    fitManager.SetRedirectOutput()  # redirect console output to files
    print(f"running {nmbProofJobs} PROOF jobs each with {nmbThreadsPerJob} threads in parallel", flush = True)
    ROOT.Proof.Go(fitManager, nmbProofJobs)
  else:
    print(f"running {nmbThreadsPerJob} threads in parallel", flush = True)
    ROOT.Here.Go(fitManager)

  fitManager.WriteThis()  # write to disk
  ROOT.gBenchmark.Show(gBenchmarkLabel)


if __name__ == "__main__":
  os.nice(18)  # run with second highest niceness level
  ROOT.gROOT.SetBatch(True)
  makePlots.setupPlotStyle()
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable
  ROOT.gBenchmark.Start("Total execution time")

  # dataset           = "030730",
  dataset           = "bggen_2017_01-ver03"
  dataFileName      = f"../ReactionEfficiency/pippippimpimpmiss_flatTree.{dataset}.root.brufit"
  outputDirName     = "./BruFitOutput"
  kinematicBinnings = [
    # ("BeamEnergy",          9,    3.0,   12.0)
    ("MissingProtonP",     10,    0.0,    3.5)
    # ("MissingProtonTheta", 13,    0.0,   65.0)
    # ("MissingProtonPhi",   10, -180.0, +180.0)
  ]

  dataSets = {
    "Total"   : None,
    "Found"   : "(TrackFound == 1)",
    "Missing" : "(TrackFound == 0)"
  }

  #TODO calculate chi^2
  # see https://root.cern/doc/master/rf109__chi2residpull_8py.html
  # and https://root-forum.cern.ch/t/how-to-correctly-extract-the-chi2-ndf-p-value-of-a-roofitresult/45956
  for dataSetName, cut in dataSets.items():
    cut = "(IsSignal == 1)" + ("" if cut is None else f" && ({cut})")
    # fit overall distribution
    performFit(
      dataFileName,
      f"{outputDirName}/{dataSetName}",
      kinematicBinnings = None,
      cut = cut,
      templateDataSigFileName = dataFileName,
      templateDataBkgFileName = dataFileName
    )
    if kinematicBinnings:
      # fit kinematic bins
      performFit(
        dataFileName,
        f"{outputDirName}/{dataSetName}",
        kinematicBinnings,
        cut,
        templateDataSigFileName = dataFileName,
        templateDataBkgFileName = dataFileName
      )

  ROOT.gBenchmark.Show("Total execution time")
