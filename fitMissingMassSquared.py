#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import os

import ROOT

import makePlots  # defines helper functions to generate histograms from data trees

makePlots.setupPlotStyle()


def defineSignalPdf(fitManager):
  # # define Gaussian signal PDF `SigPdf`
  # fitManager.SetUp().FactoryPDF(f"Gaussian::SigPdf({fitVariable}, mean_SigPdf[{meanStartVal}, 0, 2], width_SigPdf[{widthStartVal}, 0.1, 3])")

  # define double Gaussian signal PDF `SigPdf`
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


def defineBackgroundPdf(fitManager):
  # # define 2nd-order positive definite polynomial as background PDF `BgPdf`
  # fitManager.SetUp().FactoryPDF(f"GenericPdf::BgPdf('@1 * @1 + (@2 + @3 * @0) * (@2 + @3 * @0)', {{{fitVariable}, p0_BgPdf[0, -100, 100], p1_BgPdf[0, -100, 100], p2_BgPdf[0, -100, 100]}})")

  # # define 2nd-order Chebychev polynomial as background PDF `BgPdf`
  # fitManager.SetUp().FactoryPDF(f"Chebychev::BgPdf({fitVariable}, {{p0_BgPdf[0, -1, 1], p1_BgPdf[0, -1, 1], p2_BgPdf[0, -1, 1]}})")

  # define 2nd-order Bernstein polynomial as background PDF `BgPdf`
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


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable
  # ROOT.EnableImplicitMT(10)  # activate implicit multi-threading for RooFit; disable using ROOT.DisableImplicitMT()

  # dataset         = "030730"
  dataset          = "bggen_2017_01-ver03"
  dataFileName     = f"../ReactionEfficiency/pippippimpimpmiss_flatTree.{dataset}.root.brufit"
  dataTreeName     = "pippippimpimpmiss"
  fitVariable      = "MissingMassSquared_Measured"  # this is also the name of the branches in the data tree and the template-data trees for signal and background
  fitRange         = "-0.5, 4"  # [(GeV/c)^2]
  eventIDName      = "EventID"
  outputDirName    = "BruFitOutput"
  regenBinnedTrees = False

  # create the sPlot fit manager and set the output directory for fit results, plots, and weights
  fitManager = ROOT.sPlot()
  fitManager.SetUp().SetOutDir(outputDirName)
  # define fit variable and set fit range
  fitManager.SetUp().LoadVariable(f"{fitVariable}[{fitRange}]")
  # define `eventID` as event-ID variable; data tree should have a double branch with this name containing a unique event ID number
  fitManager.SetUp().SetIDBranchName(eventIDName)

  defineSignalPdf(fitManager)
  defineBackgroundPdf(fitManager)

  # # define kinematic bins
  # fitManager.Bins().LoadBinVar("BeamEnergy", 9, 3.0, 12.0)

  # load data to be fitted
  dataAreBinned = fitManager.Bins().GetBins().GetVarAxis().size() > 0
  binFileNames = binnedTreeFilesIn(outputDirName) if dataAreBinned else None
  if binFileNames is None or regenBinnedTrees:
    if dataAreBinned:
      if not regenBinnedTrees:
        print("Could not find (all) binned tree files")
      print("Regenerating binned tree files")
    fitManager.LoadData(dataTreeName, dataFileName, "Data")
  else:
    print("Using existing binned tree files:")
    for binFileName in binFileNames:
      print(f"    {binFileName}")
    fitManager.ReloadData(dataTreeName, dataFileName, "Data")

  # set fit options
  # see https://root.cern/doc/master/classRooAbsPdf.html#a52c4a5926a161bcb72eab46890b0590e
  fitManager.SetUp().AddFitOption(ROOT.RooFit.BatchMode(True))  # computes a batch of likelihood values at a time, uses faster math functions and possibly auto vectorization
                                                                # !Note! RooBatchCompute Library was revamped in ROOT 6.26/00 see https://github.com/root-project/root/tree/master/roofit/batchcompute
  fitManager.SetUp().AddFitOption(ROOT.RooFit.NumCPU(10))       # parallelizes calculation of likelihood using the given number of cores
  # fitManager.SetUp().AddFitOption(ROOT.RooFit.Parallelize(10))  # ROOT 6.28/00 global parallelization settings: enables use of RooFit's parallel minimization backend using the given number of cores
  fitManager.SetUp().AddFitOption(ROOT.RooFit.PrintLevel(2))
  fitManager.SetUp().AddFitOption(ROOT.RooFit.Timer(True))      # times CPU and wall clock consumption of fit steps
  # perform fit an plot fit result
  ROOT.Here.Go(fitManager)
  # fitManager.SetRedirectOutput()  # redirect console output to files
  # # ln -s /group/halld/Software/builds/Linux_CentOS7.7-x86_64-gcc4.8.5/root/root-6.24.04/lib/RooStats.pcm /u/home/bgrube/.proof/cache/libRooStats_rdict.pcm
  # ROOT.Proof.Go(fitManager, 9)
