#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import os

import ROOT

import makePlots  # defines helper functions to generate histograms from data trees

makePlots.setupPlotStyle()


def defineSignalPdf(
  fitManager,
  fitVariable
):
  print("Defining signal PDF 'SigPdf'", flush = True)
  # # define Gaussian
  # fitManager.SetUp().FactoryPDF(f"Gaussian::SigPdf({fitVariable}, mean_SigPdf[{meanStartVal}, 0, 2], width_SigPdf[{widthStartVal}, 0.1, 3])")

  # define double Gaussian
  # # separate means
  # fitManager.SetUp().FactoryPDF("SUM::SigPdf("
  #   f"r_SigPdf[0.5, 0, 1] * Gaussian::SigPdf_N1({fitVariable}, mean1_SigPdf[1.0,         0, 2], width1_SigPdf[1.0, 0.01, 2]),"  # wide Gaussian
  #                         f"Gaussian::SigPdf_N2({fitVariable}, mean2_SigPdf[{0.9383**2}, 0, 2], width2_SigPdf[0.2, 0.01, 2])"   # narrow Gaussian
  #   ")")
  # same mean
  fitManager.SetUp().FactoryPDF("SUM::SigPdf("
    f"r_SigPdf[0.5, 0, 1] * Gaussian::SigPdf_N1({fitVariable}, mean_SigPdf[{0.9383**2}, 0, 2], width1_SigPdf[1.0, 0.01, 2]),"  # wide Gaussian
                          f"Gaussian::SigPdf_N2({fitVariable}, mean_SigPdf,                    width2_SigPdf[0.2, 0.01, 2])"   # narrow Gaussian
    ")")

  sigPdfWeightStartVal = 1.0
  fitManager.SetUp().LoadSpeciesPDF("SigPdf", sigPdfWeightStartVal)


def defineBackgroundPdf(
  fitManager,
  fitVariable
):
  #!TODO rename Bg -> Bkg
  print("Defining background PDF 'BgPdf'", flush = True)
  # # define 2nd-order positive-definite polynomial
  # fitManager.SetUp().FactoryPDF(f"GenericPdf::BgPdf('@1 * @1 + (@2 + @3 * @0) * (@2 + @3 * @0)', {{{fitVariable}, p0_BgPdf[0, -100, 100], p1_BgPdf[0, -100, 100], p2_BgPdf[0, -100, 100]}})")

  # # define 2nd-order Chebychev polynomial
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BgPdf({fitVariable}, {{p0_BgPdf[0, -1, 1], p1_BgPdf[0, -1, 1], p2_BgPdf[0, -1, 1]}})")

  # define 2nd-order Bernstein polynomial
  # see https://root.cern.ch/doc/master/classRooBernstein.html
  fitManager.SetUp().FactoryPDF(f"Bernstein::BgPdf({fitVariable}, {{p0_BgPdf[0, 0, 1], p1_BgPdf[0, 0, 1], p2_BgPdf[0, 0, 1]}})")

  sigPdfWeightStartVal = 1.0
  fitManager.SetUp().LoadSpeciesPDF("BgPdf", sigPdfWeightStartVal)


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


def readSidebandWeights(
  inputFileName,
  inputTreeName,
  weightBranchName,
  comboIdName,
  sWeightFileName,
  sWeightLabel,  # label used to access weights
  cut = None
):
  print(f"Reading weights '{weightBranchName}' from tree '{inputTreeName}' in file '{inputFileName}'"
  f" and writing them to '{sWeightFileName}' with label '{sWeightLabel}'"
  + (f" while applying cut(s) '{cut}'" if not cut is None else ""), flush = True)
  currentDir = ROOT.gDirectory
  inputFile = ROOT.TFile.Open(inputFileName, "READ")
  inputTree = inputFile.Get(inputTreeName)
  currentDir.cd()
  weights = ROOT.Weights("HSsWeights")  # name of the Weights object
  weights.SetFile(sWeightFileName)
  weights.SetSpecies(sWeightLabel)
  weights.SetIDName(comboIdName)
  weights.WeightBySelection(inputTree, ("(1)" if cut is None else cut), weightBranchName)
  weights.SortWeights()
  weights.Save()


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


def performFit(
  dataFileName,
  outputDirName,
  kinematicBinnings,
  cut               = None,
  dataTreeName      = "pippippimpimpmiss",
  fitVariable       = "MissingMassSquared_Measured",  # this is also the name of the branches in the data tree and the template-data trees for signal and background
  fitRange          = "-0.5, 4.0",  # [(GeV/c)^2]
  comboIdName       = "ComboID",
  regenBinnedTrees  = False,
  nmbThreadsPerJob  = 5,
  nmbProofJobs      = 9
):
  print(f"fitting data in '{dataFileName}'" + ("" if cut is None else f" applying cut '{cut}'")
    + f" using binning '{kinematicBinnings}' and writing output to '{outputDirName}'")
  print(f"reading fit variable '{fitVariable}' from tree '{dataTreeName}' and using fit range {fitRange}")
  gBenchmarkLabel = f"Total time for fit in '{outputDirName}'"
  ROOT.gBenchmark.Start(gBenchmarkLabel)

  # create the fit manager and set the output directory for fit results, plots, and weights
  fitManager = ROOT.FitManager()
  fitManager.SetUp().SetOutDir(outputDirName)
  # define fit variable and set fit range
  fitManager.SetUp().LoadVariable(f"{fitVariable}[{fitRange}]")
  # define combo-ID variable
  # the data tree must have a double branch of the given name containing a unique combo-ID number
  fitManager.SetUp().SetIDBranchName(comboIdName)

  # define fit model
  defineSignalPdf(fitManager, fitVariable)
  defineBackgroundPdf(fitManager, fitVariable)

  # define kinematic bins
  for binning in kinematicBinnings:
    fitManager.Bins().LoadBinVar(*binning)

  # create RF-sideband weights for data to be fitted
  rfSWeightFileName = f"{outputDirName}/sidebandWeightsData.root"
  rfSWeightLabel    = "RfSideband"
  readSidebandWeights(
    inputFileName    = dataFileName,
    inputTreeName    = dataTreeName,
    weightBranchName = "AccidWeightFactor",
    comboIdName      = comboIdName,
    sWeightFileName  = rfSWeightFileName,
    sWeightLabel     = rfSWeightLabel,
    cut              = cut
  )
  # apply weights for RF-sideband subtraction
  fitManager.Data().LoadWeights(rfSWeightLabel, rfSWeightFileName)

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

  #TODO add constraints for fugde parameters

  # perform fit an plot fit result
  setRooFitOptions(fitManager, nmbThreadsPerJob)
  if kinematicBinnings:
    fitManager.SetRedirectOutput()  # redirect console output to files
    print(f"running {nmbProofJobs} each with {nmbThreadsPerJob} threads in parallel", flush = True)
    ROOT.Proof.Go(fitManager, nmbProofJobs)
  else:
    print(f"running {nmbThreadsPerJob} threads in parallel", flush = True)
    ROOT.Here.Go(fitManager)

  fitManager.WriteThis()  # write to disk
  ROOT.gBenchmark.Show(gBenchmarkLabel)


if __name__ == "__main__":
  os.nice(18)  # run with second highest niceness level
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable

  # dataset           = "030730",
  dataset           = "bggen_2017_01-ver03"
  dataFileName      = f"../ReactionEfficiency/pippippimpimpmiss_flatTree.{dataset}.root.brufit"
  outputDirName     = "BruFitOutput"
  kinematicBinnings = [  # list of tuples [ (variable, nmb of bins, min value, max value) ]
    ("BeamEnergy",       9,  3.0, 12.0)
    # ("MissingProtonPhi", 3, -180, +180)
  ]

  dataSets = {
    "Total"   : None,
    "Found"   : "(TrackFound == 1)",
    "Missing" : "(TrackFound == 0)"
  }

  for dataSetName, cut in dataSets.items():
    performFit(
      dataFileName,
      f"{outputDirName}/{dataSetName}",
      kinematicBinnings,
      cut
    )
