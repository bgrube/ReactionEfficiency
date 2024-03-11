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
  """Holds setup information for a fit"""
  dataFileName:     str
  templateFileName: str
  commonCut:        str
  dataCut:          str
  useSigPdf:        bool
  useBkgPdf:        bool
  fitDirName:       str


def createPdf(
  fitManager:       "ROOT.FitManager",
  fitSetup:         FitSetup,
  treeName:         str,
  fitVariable:      str,
  weightLabel:      str,
  weightBranchName: str,
  comboIdName:      str,
  pdfName:          str,
  cut:              str,
) -> None:
  """Creates a template histogram PDF from weighted Monte Carlo data"""
  # create RF-sideband weights for background histogram PDF
  weightObjectName = f"{weightLabel}Weights{pdfName}"
  weightFileName   = f"{fitSetup.fitDirName}/{weightObjectName}.root"
  print(f">>> Reading weights '{weightBranchName}' from tree '{treeName}' in file '{fitSetup.templateFileName}'"
  f" and  writing them to key '{weightObjectName}' in file '{weightFileName}', and assigning label '{weightObjectName}' while applying cut(s) '{cut}'")
  currentDir = ROOT.gDirectory
  inputFile = ROOT.TFile.Open(fitSetup.templateFileName, "READ")
  inputTree = inputFile.Get(treeName)
  currentDir.cd()
  weights = ROOT.Weights(weightObjectName)  # set name of the Weights object
  weights.SetFile(weightFileName)
  weights.SetSpecies(weightObjectName)
  weights.SetIDName(comboIdName)
  weights.WeightBySelection(inputTree, cut, weightBranchName)
  weights.SortWeights()
  weights.Save()
  # define histogram PDF with fudge parameters
  fitManager.SetUp().FactoryPDF(
    f"RooHSEventsHistPDF::{pdfName}({fitVariable}, smear_{pdfName}[0], shift_{pdfName}[0], scale_{pdfName}[1])"
    # f"RooHSEventsHistPDF::{pdfName}({fitVariable}, smear_{pdfName}[0, 0, 0], shift_{pdfName}[0, 0, 0], scale_{pdfName}[1, 1, 1])"
    # f"RooHSEventsHistPDF::{pdfName}({fitVariable}, smear_{pdfName}[0, 0, 0.5], shift_{pdfName}[0, -0.25, 0.25], scale_{pdfName}[1, 0.5, 1.5])"
    + f"WEIGHTS@{weightObjectName},{weightFileName},{weightObjectName}"  # apply sWeights created above; !Note! no whitespace allowed in this string
  )
  # constrain PDF fudge parameters, if left free
  # pdf = fitManager.SetUp().WS().pdf(pdfName)
  # fitManager.SetUp().AddGausConstraint(pdf.AlphaConstraint())  # constrain smear
  # fitManager.SetUp().AddGausConstraint(pdf.OffConstraint())    # constrain shift
  # fitManager.SetUp().AddGausConstraint(pdf.ScaleConstraint())  # constrain scale
  # load data for template histograms; ensure that calling code has already defined the kinematical binnings
  print(f">>> Loading data for histogram PDF '{pdfName}' from tree '{treeName}' in file '{fitSetup.templateFileName}'"
  + f" and applying weights from '{weightObjectName}' in file '{weightFileName}'")
  fitManager.LoadSimulated(treeName, fitSetup.templateFileName, pdfName)
  fitManager.SetUp().LoadSpeciesPDF(pdfName, 1.0)


if __name__ == "__main__":

  bggenFileName = "../data/MCbggen/2017_01-ver03_goodToF/pippippimpimpmiss_flatTree.MCbggen_2017_01-ver03_goodToF.root.brufit"
  # define fit setups
  fitSetups: list[FitSetup] = [
    # fit real data with bggen templates
    FitSetup(
      dataFileName     = "../data/RD/2017_01-ver03_goodToF/pippippimpimpmiss_flatTree.RD_2017_01-ver03_goodToF.root.brufit",
      templateFileName = bggenFileName,
      commonCut        = "(NmbUnusedShowers == 0) && (MissingProtonP > 0.5) && ((10 < MissingProtonTheta) && (MissingProtonTheta < 20)) && ((4.5 < MissingProtonP) && (MissingProtonP < 9))",
      dataCut          = "(1)",
      useSigPdf        = True,
      useBkgPdf        = True,
      fitDirName       = "./out",
    ),
    # fit bggen signal data with true template; should match perfectly
    FitSetup(
      dataFileName     = bggenFileName,
      templateFileName = bggenFileName,
      commonCut        = "(1)",
      dataCut          = "(IsSignal == 1)",  # cut applied only to data
      useSigPdf        = True,  # en/disables signal PDF
      useBkgPdf        = False,  # en/disables background PDF
      fitDirName       = "./outTruthSig",
    ),
    # fit bggen background data with true template; should match perfectly
    FitSetup(
      dataFileName     = bggenFileName,
      templateFileName = bggenFileName,
      commonCut        = "(1)",
      dataCut          = "(IsSignal == 0)",  # cut applied only to data
      useSigPdf        = False,  # en/disables signal PDF
      useBkgPdf        = True,  # en/disables background PDF
      fitDirName       = "./outTruthBkg",
    ),
    # fit bggen signal + background data with true templates; should match perfectly
    FitSetup(
      dataFileName     = bggenFileName,
      templateFileName = bggenFileName,
      commonCut        = "(1)",
      dataCut          = "(1)",  # cut applied only to data
      useSigPdf        = True,  # en/disables signal PDF
      useBkgPdf        = True,  # en/disables background PDF
      fitDirName       = "./outTruth",
    ),
  ]

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

  for fitSetup in fitSetups:
    print(f">>> fitting data from file '{fitSetup.dataFileName}' with templates from file '{fitSetup.templateFileName}' "
          f"using common cut '{fitSetup.commonCut}', data cut '{fitSetup.dataCut}', "
          f"signal PDF = {fitSetup.useSigPdf}, and background PDF = {fitSetup.useBkgPdf}; "
          f"writing fit results to directory '{fitSetup.fitDirName}'")

    sigCut = f"(({fitSetup.commonCut})) && ((IsSignal == 1))"  # cut that selects signal events
    bkgCut = f"(({fitSetup.commonCut})) && ((IsSignal == 0))"  # cut that selects background events

    # create the fit manager and set the output directory for fit results, plots, and weights
    fitManager = ROOT.FitManager()
    fitManager.SetUp().SetOutDir(fitSetup.fitDirName)

    # define fit variable and set fit range
    print(f">>> Reading fit variable '{fitVariable}' from tree '{treeName}' and using fit range {fitRange}")
    fitManager.SetUp().LoadVariable(f"{fitVariable}[{fitRange}]")
    # define combo-ID variable
    # the data tree must have a double branch of the given name containing a unique combo-ID number
    fitManager.SetUp().SetIDBranchName(comboIdName)

    # create fit function
    if fitSetup.useSigPdf:
      # signal histogram PDF
      createPdf(fitManager, fitSetup, treeName, fitVariable, weightLabel, weightBranchName, comboIdName, sigPdfName, sigCut)
    if fitSetup.useBkgPdf:
      # background PDF
      createPdf(fitManager, fitSetup, treeName, fitVariable, weightLabel, weightBranchName, comboIdName, bkgPdfName, bkgCut)

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
