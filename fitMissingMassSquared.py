#!/usr/bin/python3
# !NOTE! only on ifarm the shebang selects the correct Python3 version for ROOT


import ROOT

import makePlots  # defines helper functions to generate histograms from data trees

makePlots.setupPlotStyle()


if __name__ == "__main__":
  ROOT.gROOT.SetBatch(True)
  ROOT.gROOT.ProcessLine(".x ~/Analysis/brufit/macros/LoadBru.C")  #TODO use BRUFIT environment variable

  # dataset      = "030730"
  dataset       = "bggen_2017_01-ver03"
  dataFileName  = f"../ReactionEfficiency/pippippimpimpmiss_flatTree.{dataset}.root.brufit"
  dataTreeName  = "pippippimpimpmiss"
  fitVariable   = "MissingMassSquared_Measured"  # this is also the name of the branches in the data tree and the template-data trees for signal and background
  fitRange      = "-0.5, 4"  # [(GeV/c)^2]
  eventIDName   = "EventID"
  outputDirName = "BruFitOutput"

  # create the sPlot fit manager and set the output directory for fit results, plots, and weights
  fitManager = ROOT.sPlot()
  fitManager.SetUp().SetOutDir(outputDirName)
  # define fit variable and set fit range
  fitManager.SetUp().LoadVariable(f"{fitVariable}[{fitRange}]")
  # define `eventID` as event-ID variable; data tree should have a double branch with this name containing a unique event ID number
  fitManager.SetUp().SetIDBranchName(eventIDName)

  # meanStartVal      = 0.9383**2  # (proton mass)^2 [GeV^2]
  # widthStartVal     = 0.5        # [GeV^2]
  pdfWeightStartVal = 1
  # # define Gaussian signal PDF `SigPdf`
  # fitManager.SetUp().FactoryPDF(f"Gaussian::SigPdf({fitVariable}, mean_SigPdf[{meanStartVal}, 0, 2], width_SigPdf[{widthStartVal}, 0.1, 3])")
  # define double Gaussian signal PDF `SigPdf`
  fitManager.SetUp().FactoryPDF("SUM::SigPdf("
    f"r_SigPdf[0.57, 0, 1] * Gaussian::SigPdf_N1({fitVariable}, mean1_SigPdf[1.0, 0, 2], width1_SigPdf[0.66, 0.01, 2]),"
                           f"Gaussian::SigPdf_N2({fitVariable}, mean2_SigPdf[0.9, 0, 2], width2_SigPdf[0.24, 0.01, 2])"
    ")")
  fitManager.SetUp().LoadSpeciesPDF("SigPdf", pdfWeightStartVal)

  # define 2nd-order Bernstein polynomial as background PDF `BgPdf`
  # see https://root.cern.ch/doc/master/classRooBernstein.html
  fitManager.SetUp().FactoryPDF(f"Bernstein::BgPdf({fitVariable}, {{p0_BgPdf[0, 0, 1], p1_BgPdf[0, 0, 1], p2_BgPdf[0, 0, 1]}})")
  fitManager.SetUp().LoadSpeciesPDF("BgPdf", pdfWeightStartVal)

  # load data to be fitted
  fitManager.LoadData(dataTreeName, dataFileName)

  # perform fit an plot fit result
  ROOT.Here.Go(fitManager)
