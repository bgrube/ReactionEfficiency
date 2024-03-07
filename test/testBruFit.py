#!/usr/bin/env python3


from dataclasses import dataclass
import functools
import subprocess

import ROOT

import plotTools

# always flush print() to reduce garbling of log files due to buffering
print = functools.partial(print, flush = True)


@dataclass
class FitSetup:
  dataFileName:     str
  templateFileName: str
  commonCut:        str
  dataCut:          str
  useSigPdf:        bool
  useBkgPdf:        bool
  fitDirName:       str



if __name__ == "__main__":

  # # fit real data with bggen templates
  # fitSetup = FitSetup(
  #   dataFileName     = "../data/RD/2017_01-ver03_goodToF/pippippimpimpmiss_flatTree.RD_2017_01-ver03_goodToF.root.brufit",
  #   templateFileName = "../data/MCbggen/2017_01-ver03_goodToF/pippippimpimpmiss_flatTree.MCbggen_2017_01-ver03_goodToF.root.brufit",
  #   commonCut        = "(NmbUnusedShowers == 0) && (MissingProtonP > 0.5) && ((10 < MissingProtonTheta) && (MissingProtonTheta < 20)) && ((4.5 < MissingProtonP) && (MissingProtonP < 9))",
  #   dataCut          = "(1)",
  #   useSigPdf        = True,
  #   useBkgPdf        = True,
  #   fitDirName       = "./out",
  # )

  # fit bggen data with true template; should match perfectly
  dataFileName = "../data/MCbggen/2017_01-ver03_goodToF/pippippimpimpmiss_flatTree.MCbggen_2017_01-ver03_goodToF.root.brufit"
  fitSetup = FitSetup(
    dataFileName     = dataFileName,
    templateFileName = dataFileName,
    commonCut        = "(1)",
    dataCut          = "(IsSignal == 1)",  # cut applied only to data
    useSigPdf        = True,  # en/disables signal PDF
    useBkgPdf        = False,  # en/disables background PDF
    # dataCut          = "(IsSignal == 0)",  # cut applied only to data
    # useSigPdf        = False,  # en/disables signal PDF
    # useBkgPdf        = True,  # en/disables background PDF
    # dataCut          = "(1)",  # cut applied only to data
    # useSigPdf        = True,  # en/disables signal PDF
    # useBkgPdf        = True,  # en/disables background PDF
    fitDirName       = "./out",
  )

  sigCut               = f"(({fitSetup.commonCut})) && ((IsSignal == 1))"  # cut that selects signal events
  bkgCut               = f"(({fitSetup.commonCut})) && ((IsSignal == 0))"  # cut that selects background events
  treeName             = "pippippimpimpmiss"  # name of tree that holds the data to fit
  fitVariable          = "MissingMassSquared_Measured"  # name of branch that holds data to fit and template-data for signal and background, respectively
  fitRange             = "-0.25, 3.75"  # [(GeV/c)^2]
  comboIdName          = "ComboID"  # name of branch with unique combo ID
  sigPdfName           = "SigPdf"
  bkgPdfName           = "BkgPdf"
  weightBranchName     = "AccidWeightFactor"
  weightLabel          = "RfSideband"

  ROOT.gROOT.SetBatch(True)
  plotTools.setupPlotStyle("../rootlogon.C")
  ROOT.gROOT.ProcessLine(".x $BRUFIT/macros/LoadBru.C")

  # create the fit manager and set the output directory for fit results, plots, and weights
  fitManager = ROOT.FitManager()
  fitManager.SetUp().SetOutDir(fitSetup.fitDirName)

  # define fit variable and set fit range
  print(f">>> Reading fit variable '{fitVariable}' from tree '{treeName}' and using fit range {fitRange}")
  fitManager.SetUp().LoadVariable(f"{fitVariable}[{fitRange}]")
  # define combo-ID variable
  # the data tree must have a double branch of the given name containing a unique combo-ID number
  fitManager.SetUp().SetIDBranchName(comboIdName)

  if fitSetup.useSigPdf:
    # 1) signal histogram PDF
    # create RF-sideband weights for signal histogram PDF
    sigWeightObjectName = f"{weightLabel}Weights{sigPdfName}"
    sigWeightFileName   = f"{fitSetup.fitDirName}/{sigWeightObjectName}.root"
    print(f">>> Reading weights '{weightBranchName}' from tree '{treeName}' in file '{fitSetup.templateFileName}'"
          f" and  writing them to key '{sigWeightObjectName}' in file '{sigWeightFileName}', and assigning label '{sigWeightObjectName}' while applying cut(s) '{sigCut}'")
    currentDir = ROOT.gDirectory
    inputFile = ROOT.TFile.Open(fitSetup.templateFileName, "READ")
    inputTree = inputFile.Get(treeName)
    currentDir.cd()
    weights = ROOT.Weights(sigWeightObjectName)  # set name of the Weights object
    weights.SetFile(sigWeightFileName)
    weights.SetSpecies(sigWeightObjectName)
    weights.SetIDName(comboIdName)
    weights.WeightBySelection(inputTree, sigCut, weightBranchName)
    weights.SortWeights()
    weights.Save()
    # define signal histogram PDF with fudge parameters
    fitManager.SetUp().FactoryPDF(
      f"RooHSEventsHistPDF::{sigPdfName}({fitVariable}, smear_{sigPdfName}[0], shift_{sigPdfName}[0], scale_{sigPdfName}[1])"
      # f"RooHSEventsHistPDF::{sigPdfName}({fitVariable}, smear_{sigPdfName}[0, 0, 0], shift_{sigPdfName}[0, 0, 0], scale_{sigPdfName}[1, 1, 1])"
      # f"RooHSEventsHistPDF::{sigPdfName}({fitVariable}, smear_{sigPdfName}[0, 0, 0.5], shift_{sigPdfName}[0, -0.25, 0.25], scale_{sigPdfName}[1, 0.5, 1.5])"
      + f"WEIGHTS@{sigWeightObjectName},{sigWeightFileName},{sigWeightObjectName}"  # apply sWeights created above; !Note! no whitespace allowed in this string
    )
    # constrain signal PDF fudge parameters, if left free
    # sigPdf = fitManager.SetUp().WS().pdf(sigPdfName)
    # fitManager.SetUp().AddGausConstraint(sigPdf.AlphaConstraint())  # constrain smear
    # fitManager.SetUp().AddGausConstraint(sigPdf.OffConstraint())    # constrain shift
    # fitManager.SetUp().AddGausConstraint(sigPdf.ScaleConstraint())  # constrain scale
    # load data for template histograms; ensure that calling code has already defined the kinematical binnings
    print(f">>> Loading data for histogram PDF '{sigPdfName}' from tree '{treeName}' in file '{fitSetup.templateFileName}'"
    + f" and applying weights from '{sigWeightObjectName}' in file '{sigWeightFileName}'")
    fitManager.LoadSimulated(treeName, fitSetup.templateFileName, sigPdfName)
    fitManager.SetUp().LoadSpeciesPDF(sigPdfName, 1.0)

  if fitSetup.useBkgPdf:
    # 2) background PDF
    # create RF-sideband weights for background histogram PDF
    bkgWeightObjectName = f"{weightLabel}Weights{bkgPdfName}"
    bkgWeightFileName   = f"{fitSetup.fitDirName}/{bkgWeightObjectName}.root"
    print(f">>> Reading weights '{weightBranchName}' from tree '{treeName}' in file '{fitSetup.templateFileName}'"
          f" and  writing them to key '{bkgWeightObjectName}' in file '{bkgWeightFileName}', and assigning label '{bkgWeightObjectName}' while applying cut(s) '{bkgCut}'")
    currentDir = ROOT.gDirectory
    inputFile = ROOT.TFile.Open(fitSetup.templateFileName, "READ")
    inputTree = inputFile.Get(treeName)
    currentDir.cd()
    weights = ROOT.Weights(bkgWeightObjectName)  # set name of the Weights object
    weights.SetFile(bkgWeightFileName)
    weights.SetSpecies(bkgWeightObjectName)
    weights.SetIDName(comboIdName)
    weights.WeightBySelection(inputTree, bkgCut, weightBranchName)
    weights.SortWeights()
    weights.Save()
    # define background histogram PDF with fudge parameters
    fitManager.SetUp().FactoryPDF(
      f"RooHSEventsHistPDF::{bkgPdfName}({fitVariable}, smear_{bkgPdfName}[0], shift_{bkgPdfName}[0], scale_{bkgPdfName}[1])"
      # f"RooHSEventsHistPDF::{bkgPdfName}({fitVariable}, smear_{bkgPdfName}[0, 0, 0.5], shift_{bkgPdfName}[0, -0.25, 0.25], scale_{bkgPdfName}[1, 0.5, 1.5])"
      + f"WEIGHTS@{bkgWeightObjectName},{bkgWeightFileName},{bkgWeightObjectName}"  # apply sWeights created above; !Note! no whitespace allowed in this string
    )
    # constrain background PDF fudge parameters, if left free
    # bkgPdf = fitManager.SetUp().WS().pdf(bkgPdfName)
    # fitManager.SetUp().AddGausConstraint(bkgPdf.AlphaConstraint())  # constrain smear
    # fitManager.SetUp().AddGausConstraint(bkgPdf.OffConstraint())    # constrain shift
    # fitManager.SetUp().AddGausConstraint(bkgPdf.ScaleConstraint())  # constrain scale
    # load data for template histograms; ensure that calling code has already defined the kinematical binnings
    print(f">>> Loading data for histogram PDF '{bkgPdfName}' from tree '{treeName}' in file '{fitSetup.templateFileName}'"
    + f" and applying weights from '{bkgWeightObjectName}' in file '{bkgWeightFileName}'")
    fitManager.LoadSimulated(treeName, fitSetup.templateFileName, bkgPdfName)
    fitManager.SetUp().LoadSpeciesPDF(bkgPdfName, 1.0)

  # create RF-sideband weights for data to be fitted
  dataWeightObjectName = f"{weightLabel}WeightsData"
  dataWeightFileName   = f"{fitSetup.fitDirName}/{dataWeightObjectName}.root"
  cut = f"(({fitSetup.commonCut}) && ({fitSetup.dataCut}))"
  print(f">>> Reading weights '{weightBranchName}' from tree '{treeName}' in file '{fitSetup.dataFileName}'"
        f" and  writing them to key '{dataWeightObjectName}' in file '{dataWeightFileName}', and assigning label '{dataWeightObjectName}' while applying cut(s) '{cut}'")
  currentDir = ROOT.gDirectory
  inputFile = ROOT.TFile.Open(fitSetup.dataFileName, "READ")
  inputTree = inputFile.Get(treeName)
  currentDir.cd()
  weights = ROOT.Weights(dataWeightObjectName)  # set name of the Weights object
  weights.SetFile(dataWeightFileName)
  weights.SetSpecies(dataWeightObjectName)
  weights.SetIDName(comboIdName)
  weights.WeightBySelection(inputTree, cut, weightBranchName)
  weights.SortWeights()
  weights.Save()
  # apply weights for RF-sideband subtraction
  fitManager.Data().LoadWeights(dataWeightObjectName, dataWeightFileName, dataWeightObjectName)

  # load and bin data to be fitted
  print(f">>> reading data to fit from tree '{treeName}' in file '{fitSetup.dataFileName}'")
  fitManager.LoadData(treeName, fitSetup.dataFileName)

  # perform fit and create fit-result plots
  print(">>> performing fit")
  fitManager.SetUp().AddFitOption(ROOT.RooFit.BatchMode(True))  # computes a batch of likelihood values at a time, uses faster math functions and possibly auto vectorization
  fitManager.SetUp().AddFitOption(ROOT.RooFit.NumCPU(25))  # parallelizes calculation of likelihood using the given number of cores
  fitManager.SetUp().AddFitOption(ROOT.RooFit.PrintLevel(2))
  fitManager.SetUp().AddFitOption(ROOT.RooFit.Timer(True))  # times CPU and wall clock consumption of fit steps
  # try MCMC algorithm with
  # mcmc = ROOT.BruMcmcCovariance(200, 100, 0.1, 0.23, 0.16, 0.3)
  # mcmc.TurnOffCovariance()  # BruMcmcCovariance only, do not proceed with covariance-based sampling, just perform basic stepping
  # fitManager.SetMinimiser(mcmc)
  ROOT.Here.Go(fitManager)

  print(">>> writing fit result")
  fitManager.WriteThis()

  print(">>> plotting fit result")
  # perform plotting in separate process to preserve plot style
  subprocess.run(f"./testBruFitMakePlots.py \"{fitSetup.fitDirName}\"  \"{fitVariable}\"", shell = True)
